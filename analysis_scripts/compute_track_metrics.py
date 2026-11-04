import sys
import os
import pandas as pd
import pyBigWig
import numpy as np

def get_extra_field_index(bw, field_name):
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
                if line == ')': break
                if in_body and line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        fields.append(parts[1].rstrip(';'))
            extra_fields = fields[3:]
            if field_name in extra_fields: return extra_fields.index(field_name)
    except Exception: pass
    return None

bed_file = sys.argv[1]
track_dir = sys.argv[2]
out_file = sys.argv[3]
output_track_name = sys.argv[4]

bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'gene_id', 'track_name'])
results = []

# Cache for open pyBigWig handles
bw_handles = {}

def get_bw(t_name):
    if t_name not in bw_handles:
        bw_path = f"{track_dir}/{t_name}.bw"
        bb_path = f"{track_dir}/{t_name}.bb"
        t_path = bw_path if os.path.exists(bw_path) else bb_path
        try:
            bw_handles[t_name] = (pyBigWig.open(t_path), t_path.endswith('.bb'))
        except Exception:
            bw_handles[t_name] = (None, False)
    return bw_handles[t_name]

# Determine global type (BigWig vs BigBed) based on the tracks used in this file
is_bigbed_global = False
for t in bed_df['track_name'].unique():
    _, is_bb = get_bw(t)
    if is_bb:
        is_bigbed_global = True
        break

for gene_id, group in bed_df.groupby('gene_id'):
    if not is_bigbed_global:
        # CONTINUOUS: Calculate Mean and Median across the WHOLE gene (all exons)
        gene_vals = []
        for _, row in group.iterrows():
            bw, is_bb = get_bw(row['track_name'])
            if bw is None: continue
            
            if row['chrom'] in bw.chroms():
                try:
                    vals = bw.values(row['chrom'], row['start'], row['end'])
                    gene_vals.extend([v for v in vals if not np.isnan(v)])
                except RuntimeError:
                    pass
        
        if gene_vals:
            results.append((gene_id, np.mean(gene_vals), np.median(gene_vals)))
        else:
            results.append((gene_id, np.nan, np.nan))
            
    else:
        # DISCRETE
        total_mutations = 0
        total_length = 0
        max_sample_count = 0
        weighted_count = 0
        
        for _, row in group.iterrows():
            bw, is_bb = get_bw(row['track_name'])
            if bw is None: continue
            
            total_length += (row['end'] - row['start'])
            sc_idx = get_extra_field_index(bw, 'sampleCount')
            
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
                    
        if total_length > 0:
            density_per_kb = (total_mutations / total_length) * 1000
            results.append((gene_id, total_mutations, density_per_kb, max_sample_count, weighted_count))
        else:
            results.append((gene_id, np.nan, np.nan, np.nan, np.nan))

# Close handles safely
for bw, _ in bw_handles.values():
    if bw is not None: bw.close()

# Cleanup naming for output (e.g. "PhyloCSF_Combined" -> "PhyloCSF")
final_prefix = output_track_name.replace("_Combined", "")

if not is_bigbed_global:
    res_df = pd.DataFrame(results, columns=['Geneid', f'{final_prefix}_Mean', f'{final_prefix}_Median'])
else:
    res_df = pd.DataFrame(results, columns=[
        'Geneid', f'{final_prefix}_Count', f'{final_prefix}_Density_per_kb',
        f'{final_prefix}_Max_sampleCount', f'{final_prefix}_Weighted_sampleCount'
    ])

res_df.to_csv(out_file, sep='\t', index=False)