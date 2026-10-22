#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm, pearsonr
from multiprocessing import Pool, cpu_count
from sklearn.decomposition import PCA
import warnings
import logging
import sys

# Silence warnings
logging.getLogger("sklearn").setLevel(logging.ERROR)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

def evaluate_uorf(args):
    uorf_id, cds_id, u_vals, c_vals, nv_nt, nv_t = args
    
    valid_mask = (u_vals > 0) & (c_vals > 0)
    if valid_mask.sum() < 10:
        return None
        
    u_log = np.log2(u_vals[valid_mask])
    c_log = np.log2(c_vals[valid_mask])
    
    if np.var(c_log) < 1e-8 or np.var(u_log) < 1e-8:
        r_val, p_corr = np.nan, np.nan
    else:
        r_val, p_corr = pearsonr(c_log, u_log)
        
    ratio = u_log - c_log
    
    # --- PCA Approach ---
    try:
        # Explicitly mean-center the coordinates
        u_log_c = u_log - np.mean(u_log)
        c_log_c = c_log - np.mean(c_log)
        
        # Stack into [n_samples, 2] matrix: Column 0 = uORF, Column 1 = CDS
        data_matrix = np.column_stack([u_log_c.values, c_log_c.values])
        
        # Fit PCA
        pca = PCA(n_components=2)
        pca.fit(data_matrix)
        
        # Extract the eigenvector for PC1
        v_x, v_y = pca.components_[0]
        
        # Calculate angle in radians using arctan2
        angle_rad = np.arctan2(v_y, v_x)
        
        # Convert to degrees and wrap to [0, 180) to handle bidirectional lines
        pc1_angle = np.degrees(angle_rad) % 180.0
        
        # Explained variance ratio of PC1 (0 to 1)
        pc1_evr = pca.explained_variance_ratio_[0]
        
    except Exception as e:
        print(f"[Worker Error] PCA failed for {uorf_id}: {e}", file=sys.stderr)
        return None
    
    # --- Tumor Shift Setup ---
    ratio_nt = ratio[ratio.index.isin(nv_nt)].dropna()
    ratio_t = ratio[ratio.index.isin(nv_t)].dropna()
    
    if len(ratio_nt) < 3 or len(ratio_t) < 3:
        return None
        
    med_diff = ratio_t.median() - ratio_nt.median()
    stat, pval_mw = mannwhitneyu(ratio_nt, ratio_t, alternative='two-sided')
    
    safe_pval = max(pval_mw, 1e-300)
    z_mw = norm.isf(safe_pval / 2) * np.sign(med_diff)
    
    return {
        'uorf_id': uorf_id,
        'cds_id': cds_id,
        'Pearson_r': r_val,
        'Pearson_pval': p_corr,
        'PC1_Angle': pc1_angle,
        'PC1_EVR': pc1_evr,
        'Ratio_Median_Delta': med_diff,
        'MWU_pval': pval_mw,
        'MWU_Z_score': z_mw
    }

def main():
    parser = argparse.ArgumentParser(description="PCA Angle Screen for uORFs")
    parser.add_argument('--input', required=True, help="Path to CDS_uORF_dORF_TPM_values.tsv")
    parser.add_argument('--meta', required=True, help="Path to Ribo-seq Samples metadata TSV")
    parser.add_argument('--output', required=True, help="Output TSV path")
    parser.add_argument('--chunk_id', type=int, default=0, help="0-indexed chunk ID for parallelization")
    parser.add_argument('--total_chunks', type=int, default=1, help="Total number of chunks")
    parser.add_argument('--cores', type=int, default=cpu_count(), help="CPU cores to use")
    parser.add_argument('--element', type=str, required=True, choices=['uORF', 'dORF', 'RT'], help="Element type to screen")
    args = parser.parse_args()

    # --- Metadata Parsing ---
    meta_df = pd.read_csv(args.meta, sep="\t")
    etiology_classification_dict = {
        'HCV':'viral', 'HBV':'viral', 'NASH':'non-viral', 'ASH':'non-viral', 'ASH/NASH':'non-viral'
    }
    meta_df['etiology'] = meta_df['Liver disease background'].map(etiology_classification_dict)
    meta_df['sample_id'] = meta_df['Biobank ID'].str.replace('-', '_')
    meta_df['sample_type'] = meta_df.apply(lambda x: 'NT' if '_NT' in str(x['sample_id']) else 'T', axis=1)
    
    # --- Data Parsing ---
    df = pd.read_csv(args.input, sep="\t")
    cols = list(df.columns)
    
    # --- Sample Filtering ---
    NT_cols_viral = [e for e in meta_df.loc[(meta_df['etiology']=='viral') & (meta_df['sample_type']=='NT')]['sample_id'] if e in cols]
    NT_cols_NONviral = [e for e in meta_df.loc[(meta_df['etiology']=='non-viral') & (meta_df['sample_type']=='NT')]['sample_id'] if e in cols] + ['P275_NT_35mins_normal']
    T_cols_viral = [e for e in meta_df.loc[(meta_df['etiology']=='viral') & (meta_df['sample_type']=='T')]['sample_id'] if e in cols]
    T_cols_NONviral = [e for e in meta_df.loc[(meta_df['etiology']=='non-viral') & (meta_df['sample_type']=='T')]['sample_id'] if e in cols] + ['AG043_02']
    
    all_samples = NT_cols_NONviral + NT_cols_viral + T_cols_NONviral + T_cols_viral
    nv_nt_samples, nv_t_samples = set(NT_cols_NONviral), set(T_cols_NONviral)

    print(f"Loaded {len(all_samples)} total samples for analysis.")

    # --- DataFrame Prep ---
    element_df = df[df['Element'] == args.element].copy()
    cds_df = df[df['Element'] == 'CDS'].copy()
    
    # Extract the base CDS ID assuming format like ENST00000649427_249_429
    element_df['mapped_cds_id'] = element_df['feature_id'].apply(lambda x: str(x).split('_')[0])
    cds_df = cds_df.drop_duplicates(subset=['feature_id']).set_index('feature_id')
    
    # --- Array Job Chunking ---
    chunk_indices = np.array_split(element_df.index, args.total_chunks)[args.chunk_id]
    element_df_chunk = element_df.loc[chunk_indices]

    worker_args = []
    for _, row in element_df_chunk.iterrows():
        norf_id = row['feature_id']
        cds_id = row['mapped_cds_id']
        if cds_id in cds_df.index:
            valid_samples = [s for s in all_samples if s in row.index and s in cds_df.columns]
            u_vals = row[valid_samples].astype(float)
            c_raw = cds_df.loc[cds_id, valid_samples]
            if isinstance(c_raw, pd.DataFrame): 
                c_raw = c_raw.iloc[0]
            worker_args.append((norf_id, cds_id, u_vals, c_raw.astype(float), nv_nt_samples, nv_t_samples))
            
    print(f"Chunk {args.chunk_id + 1}/{args.total_chunks}: Screening {len(worker_args)} pairs using {args.cores} cores...")
    
    # --- Parallel Execution ---
    with Pool(processes=args.cores) as pool:
        results = pool.map(evaluate_uorf, worker_args)

    res_df = pd.DataFrame([r for r in results if r is not None])
    res_df.to_csv(args.output, sep='\t', index=False)
    print(f"Saved results to {args.output}")

if __name__ == '__main__':
    main()