import sys
import os
import pandas as pd
import pyBigWig
import numpy as np

def get_extra_field_index(bw, field_name):
    """Return 0-based index of field_name in BigBed extra (non-BED3) fields via AutoSQL."""
    try:
        sql = bw.SQL()
        if sql:
            fields = []
            in_body = False
            for line in sql.split('\n'):
                line = line.strip()
                if line == '(':
                    in_body = True
                    continue
                if line == ')':
                    break
                if in_body and line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        fields.append(parts[1].rstrip(';'))
            extra_fields = fields[3:]  # skip chrom, chromStart, chromEnd
            if field_name in extra_fields:
                return extra_fields.index(field_name)
    except Exception:
        pass
    return None


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
        results.append((gene_id, np.nan, np.nan) if not is_bigbed else (gene_id, np.nan, np.nan, np.nan, np.nan))
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
        # DISCRETE: Count, Density, Max sampleCount, Weighted sampleCount for BigBeds (gnomAD, TCGA)
        total_mutations = 0
        total_length = 0
        max_sample_count = 0
        weighted_count = 0
        sc_idx = get_extra_field_index(bw, 'sampleCount') # sampleCount is present in TCGA .bb storing recurrency of mutations
        for _, row in group.iterrows():
            total_length += (row['end'] - row['start'])
            if row['chrom'] in bw.chroms():
                try:
                    entries = bw.entries(row['chrom'], row['start'], row['end'])
                    if entries is not None:
                        total_mutations += len(entries)
                        if sc_idx is not None:
                            for entry in entries:
                                try:
                                    sc = int(entry[2].split('\t')[sc_idx])
                                    max_sample_count = max(max_sample_count, sc)
                                    weighted_count += sc
                                except (ValueError, IndexError):
                                    pass
                except RuntimeError:
                    pass
        bw.close()
        
        if total_length > 0:
            density_per_kb = (total_mutations / total_length) * 1000
            results.append((gene_id, total_mutations, density_per_kb, max_sample_count, weighted_count))
        else:
            results.append((gene_id, np.nan, np.nan, np.nan, np.nan))

# Build output DataFrame with type-appropriate columns
if not is_bigbed:
    res_df = pd.DataFrame(results, columns=[
        'Geneid', f'{track_name}_Mean', f'{track_name}_Median'
    ])
else:
    res_df = pd.DataFrame(results, columns=[
        'Geneid', f'{track_name}_Count', f'{track_name}_Density_per_kb',
        f'{track_name}_Max_sampleCount', f'{track_name}_Weighted_sampleCount'
    ])

# Save output
res_df.to_csv(out_file, sep='\t', index=False)