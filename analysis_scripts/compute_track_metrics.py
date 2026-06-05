import sys
import os
import pandas as pd
import pyBigWig
import numpy as np

# Parse command line arguments
bed_file = sys.argv[1]
track_dir = sys.argv[2]
out_file = sys.argv[3]
track_name = sys.argv[4]

# Determine if track is BigWig or BigBed based on available files
bw_path = f"{track_dir}/{track_name}.bw"
bb_path = f"{track_dir}/{track_name}.bb"

track_path = bw_path if os.path.exists(bw_path) else bb_path
is_bigbed = track_path.endswith('.bb')

# Read BED file 
bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene_id', 'track_name'])
results = []

# Group by gene to pull all CDS arrays across the whole gene
for gene_id, group in bed_df.groupby('gene_id'):
    try:
        bw = pyBigWig.open(track_path)
    except Exception:
        results.append((gene_id, np.nan, np.nan))
        continue
    
    if not is_bigbed:
        # CONTINUOUS: Calculate Mean and Median for BigWigs (PhyloP, PhyloCSF)
        gene_vals = []
        for _, row in group.iterrows():
            if row['chrom'] in bw.chroms():
                try:
                    vals = bw.values(row['chrom'], row['start'], row['end'])
                    gene_vals.extend([v for v in vals if not np.isnan(v)])
                except RuntimeError:
                    pass
        bw.close()
        
        if gene_vals:
            results.append((gene_id, np.mean(gene_vals), np.median(gene_vals)))
        else:
            results.append((gene_id, np.nan, np.nan))
            
    else:
        # DISCRETE: Calculate Count and Density for BigBeds (gnomAD, TCGA)
        total_mutations = 0
        total_length = 0
        for _, row in group.iterrows():
            total_length += (row['end'] - row['start'])
            if row['chrom'] in bw.chroms():
                try:
                    entries = bw.entries(row['chrom'], row['start'], row['end'])
                    if entries is not None:
                        total_mutations += len(entries)
                except RuntimeError:
                    pass
        bw.close()
        
        if total_length > 0:
            density_per_kb = (total_mutations / total_length) * 1000
            results.append((gene_id, total_mutations, density_per_kb))
        else:
            results.append((gene_id, np.nan, np.nan))

# Dynamic column naming based on file type
col1 = 'Mean' if not is_bigbed else 'Count'
col2 = 'Median' if not is_bigbed else 'Density_per_kb'

res_df = pd.DataFrame(results, columns=['Geneid', f'{track_name}_{col1}', f'{track_name}_{col2}'])

# Save output
res_df.to_csv(out_file, sep='\t', index=False)