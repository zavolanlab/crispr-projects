#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Author : Aleksei Mironov
# Company: Mihaela Zavolan, Biozentrum, Basel
# -----------------------------------------------------------------------------

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd

def main():
    """ Get genomic segmentation out of the .gtf file"""

    __doc__ = "Get genomic segmentation out of the .gtf file"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_gtf",
                        dest="input_gtf_file_path",
                        help="Path to the input gtf annotation file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--input_genome_fai",
                        dest="input_genome_fai_file_path",
                        help="Path to the input genome fai index file",
                        required=True,
                        metavar="FILE",)
    parser.add_argument("--output_exon_intron_gs",
                        dest="output_exon_intron_gs_file_path",
                        help="path for the output file exon_intron_genome_segmentation.bed",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--output_binned_gs",
                        dest="output_binned_gs_file_path",
                        help="path for the output file binned_genome_segmentation.bed",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--bin_size",
                    dest="bin_size",
                    help="Separate segments into sub-segments of this length, related to the previous argument",
                    required=True,
                    metavar="FILE")    
    parser.add_argument("--temp_dir",
                    dest="temp_dir_path",
                    help="path to the temporary directory",
                    required=True,
                    metavar="FILE")
    parser.add_argument("--output_modified_gtf",
                        dest="output_modified_gtf_file_path",
                        help="path for the output file modified_annotation.gtf",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--gene_flank",
                        dest="gene_flank",
                        help="gene_flank to define intergenic regions flanking the genes, nt",
                        required=True,
                        metavar="FILE")    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    temp_dir_path = options.temp_dir_path+'/'
    gene_flank = int(options.gene_flank)
    bin_size = int(options.bin_size)
    
    # create output directories if they don't exist
    for file_path in [options.output_exon_intron_gs_file_path, temp_dir_path]:
        outdir = '/'.join((file_path.split('/')[:-1]))+'/'
        out = subprocess.check_output('mkdir -p '+outdir, shell=True)
    
    input_gtf = pd.read_csv(options.input_gtf_file_path,delimiter="\t",index_col=None,header=None)
    
    # rename gene_ids so that they are not too long
    input_gtf['gene_id'] = input_gtf[8].str.split('gene_id "',expand=True)[1].str.split('";',expand=True)[0]
    gene_id_df = input_gtf[['gene_id']].drop_duplicates().reset_index(drop=True)
    gene_id_df['new_gene_id'] = gene_id_df.index
    input_gtf = pd.merge(input_gtf,gene_id_df,how='inner',on='gene_id')
    input_gtf[8] = input_gtf[8].str.split('gene_id "',expand=True)[0]+'gene_id "'+\
    input_gtf['new_gene_id'].astype('str')+'";'+input_gtf[8].str.split('gene_id "',expand=True)[1].str.split('";',expand=True,n=1)[1]
    
    input_gtf = input_gtf[list(range(0,9,1))]
    input_gtf.to_csv(options.output_modified_gtf_file_path, sep=str('\t'),header=False,index=None)
    
    for sel_strand in ['+','-']:
    
        genome_fai = pd.read_csv(options.input_genome_fai_file_path,delimiter="\t",index_col=None,header=None)
        genome_fai['start'] = 0
        genome_fai['score'] = 0
        genome_fai['chr_name'] = genome_fai[0]
        genome_fai = genome_fai.loc[genome_fai['chr_name'].isin(list(input_gtf[0].unique()))].reset_index(drop=True)
        genome_fai['strand'] = sel_strand
        genome_fai[[0,'start',1,'chr_name','score','strand']].to_csv(temp_dir_path+'chromosomes.bed', sep=str('\t'),header=False,index=None)
        
        cur_input_gtf = input_gtf.loc[input_gtf[6]==sel_strand].reset_index(drop=True).copy()
        
        genes = cur_input_gtf.loc[cur_input_gtf[2]=='gene'].reset_index(drop=True)
        genes['gene_id'] = genes[8].str.split('gene_id "',expand=True)[1].str.split('";',expand=True)[0]
        genes['score'] = 0
        genes = genes[[0,3,4,'gene_id','score',6]]
        genes[3] = genes[3]-1
        genes = genes.sort_values([0,3,4],ascending=True).reset_index(drop=True)
        genes.to_csv(temp_dir_path+'genes_original.bed', sep=str('\t'),header=False,index=None)
        
        # add flanks
        genes_with_flanks = genes.copy()
        genes_with_flanks.columns = list(range(0,len(genes_with_flanks.columns)))
        genes_with_flanks['upstream_start'] = genes_with_flanks[1]-gene_flank
        genes_with_flanks['tmp'] = 0
        genes_with_flanks['upstream_start'] = genes_with_flanks[['upstream_start','tmp']].max(1)
        genes_with_flanks['upstream_end'] = genes_with_flanks[1]
        
        genes_with_flanks = pd.merge(genes_with_flanks,genome_fai[[0,1]].rename(columns={1:'chr_end'}),how='inner',on=[0])
        genes_with_flanks['downstream_start'] = genes_with_flanks[2]
        genes_with_flanks['downstream_end'] = genes_with_flanks[2]+gene_flank
        genes_with_flanks['downstream_end'] = genes_with_flanks[['downstream_end','chr_end']].min(1)
        genes_with_flanks['strand'] = sel_strand
        genes_with_flanks = genes_with_flanks.reset_index(drop=True)
        
        if sel_strand == '+':
            five_prime_gene_flanks = genes_with_flanks[[0,'upstream_start','upstream_end',3,'tmp','strand']].copy()
            three_prime_gene_flanks = genes_with_flanks[[0,'downstream_start','downstream_end',3,'tmp','strand']].copy()
        else:
            five_prime_gene_flanks = genes_with_flanks[[0,'downstream_start','downstream_end',3,'tmp','strand']].copy()
            three_prime_gene_flanks = genes_with_flanks[[0,'upstream_start','upstream_end',3,'tmp','strand']].copy()
        
        five_prime_gene_flanks.columns = list(range(0,len(five_prime_gene_flanks.columns)))
        three_prime_gene_flanks.columns = list(range(0,len(three_prime_gene_flanks.columns)))
        
        five_prime_gene_flanks = five_prime_gene_flanks.loc[five_prime_gene_flanks[1]<five_prime_gene_flanks[2]]
        three_prime_gene_flanks = three_prime_gene_flanks.loc[three_prime_gene_flanks[1]<three_prime_gene_flanks[2]]
        
        five_prime_gene_flanks.sort_values([0,1,2],ascending=True).reset_index(drop=True).to_csv(temp_dir_path+'five_prime_gene_flanks.bed', sep=str('\t'),header=False,index=None)
        three_prime_gene_flanks.sort_values([0,1,2],ascending=True).reset_index(drop=True).to_csv(temp_dir_path+'three_prime_gene_flanks.bed', sep=str('\t'),header=False,index=None)
        
        genes_with_flanks.sort_values([0,'upstream_start','downstream_end'],ascending=True).reset_index(drop=True)[[0,'upstream_start','downstream_end',3,'tmp','strand']].to_csv(temp_dir_path+'genes_with_flanks.bed', sep=str('\t'),header=False,index=None)
            
        # putting overlapping genes into disjoint regions

        # bedops requires files to be sorted with sort-bed
        command = 'sort-bed '+temp_dir_path+'genes_with_flanks.bed > '+\
        temp_dir_path+'genes_with_flanks.bed.sorted && rm '+temp_dir_path+'genes_with_flanks.bed && mv '+temp_dir_path+'genes_with_flanks.bed.sorted '+temp_dir_path+'genes_with_flanks.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedops --partition '+temp_dir_path+'genes_with_flanks.bed > '+temp_dir_path+'genes_partion.non_sorted.bed'
        out = subprocess.check_output(command, shell=True)

        command = 'bedtools sort -i '+temp_dir_path+'genes_partion.non_sorted.bed > '+temp_dir_path+'genes_partion.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedtools intersect -loj -sorted -a '+temp_dir_path+'genes_partion.bed -b '+temp_dir_path+\
        """genes_with_flanks.bed | bedtools groupby -g 1,2,3,8,9 -c 7 -o distinct | awk '{ print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}' > """+temp_dir_path+'genes_disjoint.bed'
        out = subprocess.check_output(command, shell=True)
        
        # get exons
        exons = cur_input_gtf.loc[cur_input_gtf[2]=='exon'].reset_index(drop=True)
        exons['exon_id'] = exons[8].str.split('exon_id "',expand=True)[1].str.split('";',expand=True)[0]
        exons['score'] = 0
        exons = exons[[0,3,4,'exon_id','score',6]]
        exons[3] = exons[3]-1
        exons = exons.sort_values([0,3,4],ascending=True).reset_index(drop=True)
        exons.to_csv(temp_dir_path+'exons.bed', sep=str('\t'),header=False,index=None)
        
        # get introns
        exons = cur_input_gtf.loc[cur_input_gtf[2]=='exon'].reset_index(drop=True)
        exons['transcript_id'] = exons[8].str.split('transcript_id "',expand=True)[1].str.split('";',expand=True)[0]
        exons = exons[[0,3,4,'transcript_id',6]]
        exons = exons.sort_values([0,6,'transcript_id',3,4],ascending=True).reset_index(drop=True)
        exons = exons.values
        introns = []
        intron_start, intron_end, prev_transcript = None, None, ''
        
        for exon in exons:
            if (intron_start is None) or (prev_transcript!=exon[3]):
                intron_start = exon[2]
                prev_transcript = exon[3]
            else:
                intron_end = exon[1]
                introns.append([exon[0],intron_start,intron_end,exon[3],0,exon[4]])
                intron_start = exon[2]
                prev_transcript = exon[3]
        introns = pd.DataFrame(introns)
        introns[2] = introns[2]-1
        introns = introns.sort_values([0,1,2],ascending=True).reset_index(drop=True)
        introns.to_csv(temp_dir_path+'introns.bed', sep=str('\t'),header=False,index=None)
        
        command = 'bedtools merge -i '+temp_dir_path+'exons.bed > '+temp_dir_path+'exons_collapsed.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedtools merge -i '+temp_dir_path+'introns.bed > '+temp_dir_path+'introns_collapsed.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedtools merge -i '+temp_dir_path+'five_prime_gene_flanks.bed > '+temp_dir_path+'five_prime_gene_flanks_collapsed.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedtools merge -i '+temp_dir_path+'three_prime_gene_flanks.bed > '+temp_dir_path+'three_prime_gene_flanks_collapsed.bed'
        out = subprocess.check_output(command, shell=True)
    
        command = 'bedtools intersect -a '+temp_dir_path+'five_prime_gene_flanks_collapsed.bed -b '+temp_dir_path+'three_prime_gene_flanks_collapsed.bed > '+temp_dir_path+'five_or_three_prime_gene_flanks_collapsed.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools subtract -a '+temp_dir_path+'five_prime_gene_flanks_collapsed.bed -b '+temp_dir_path+'five_or_three_prime_gene_flanks_collapsed.bed > '+temp_dir_path+'five_prime_only_gene_flanks_collapsed.bed'
        out = subprocess.check_output(command, shell=True)    
        command = 'bedtools subtract -a '+temp_dir_path+'three_prime_gene_flanks_collapsed.bed -b '+temp_dir_path+'five_or_three_prime_gene_flanks_collapsed.bed > '+temp_dir_path+'three_prime_only_gene_flanks_collapsed.bed'
        out = subprocess.check_output(command, shell=True)    
        
        
        command = 'bedtools subtract -a '+temp_dir_path+'genes_disjoint.bed -b '+temp_dir_path+'introns_collapsed.bed > '+temp_dir_path+'temp.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools subtract -a '+temp_dir_path+'temp.bed -b '+temp_dir_path+'exons_collapsed.bed > '+temp_dir_path+'gene_flanks.bed'
        out = subprocess.check_output(command, shell=True)
    
        command = 'bedtools intersect -a '+temp_dir_path+'gene_flanks.bed -b '+temp_dir_path+'five_or_three_prime_gene_flanks_collapsed.bed | '+\
        """bedtools groupby -g 1,2,3,5,6 -c 4 -o distinct | awk '{ print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}' > """+temp_dir_path+'five_or_three_prime_gene_flanks_final.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools intersect -a '+temp_dir_path+'gene_flanks.bed -b '+temp_dir_path+'five_prime_only_gene_flanks_collapsed.bed | '+\
        """bedtools groupby -g 1,2,3,5,6 -c 4 -o distinct | awk '{ print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}' > """+temp_dir_path+'five_prime_only_gene_flanks_final.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools intersect -a '+temp_dir_path+'gene_flanks.bed -b '+temp_dir_path+'three_prime_only_gene_flanks_collapsed.bed | '+\
        """bedtools groupby -g 1,2,3,5,6 -c 4 -o distinct | awk '{ print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}' > """+temp_dir_path+'three_prime_only_gene_flanks_final.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedtools merge -i '+temp_dir_path+'gene_flanks.bed > '+temp_dir_path+'gene_flanks_collapsed.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools subtract -a '+temp_dir_path+'genes_disjoint.bed -b '+temp_dir_path+'gene_flanks_collapsed.bed > '+temp_dir_path+'exonic_and_intronic_parts.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedtools subtract -a '+temp_dir_path+'exonic_and_intronic_parts.bed -b '+temp_dir_path+'introns_collapsed.bed > '+temp_dir_path+'alw_exonic.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedtools subtract -a '+temp_dir_path+'exonic_and_intronic_parts.bed -b '+temp_dir_path+'exons_collapsed.bed > '+temp_dir_path+'alw_intronic.bed'
        out = subprocess.check_output(command, shell=True)
        
        command = 'bedtools subtract -a '+temp_dir_path+'exonic_and_intronic_parts.bed -b '+temp_dir_path+'alw_exonic.bed > '+temp_dir_path+'temp.bed'
        out = subprocess.check_output(command, shell=True)
        command = 'bedtools subtract -a '+temp_dir_path+'temp.bed -b '+temp_dir_path+'alw_intronic.bed > '+temp_dir_path+'exonic_or_intronic.bed'
        out = subprocess.check_output(command, shell=True)
    
        non_intergenic_files = ['exonic_or_intronic','alw_intronic','alw_exonic','five_or_three_prime_gene_flanks_final','five_prime_only_gene_flanks_final','three_prime_only_gene_flanks_final']
        non_intergenic_files = list(temp_dir_path+pd.Series(non_intergenic_files)+'.bed')
        
        command = 'cat '+' '.join(non_intergenic_files)+' | bedtools sort > '+temp_dir_path+'non_intergenic.bed'
        out = subprocess.check_output(command, shell=True)
    
        command = 'bedtools merge -i '+temp_dir_path+'non_intergenic.bed > '+temp_dir_path+'non_intergenic_collapsed.bed'
        out = subprocess.check_output(command, shell=True)    
        
        command = 'bedtools subtract -a '+temp_dir_path+'chromosomes.bed -b '+temp_dir_path+'non_intergenic_collapsed.bed > '+temp_dir_path+'intergenic.bed'
        out = subprocess.check_output(command, shell=True)
    
        exonic_or_intronic = pd.read_csv(temp_dir_path+'exonic_or_intronic.bed',delimiter="\t",index_col=None,header=None)
        alw_exonic = pd.read_csv(temp_dir_path+'alw_exonic.bed',delimiter="\t",index_col=None,header=None)
        alw_intronic = pd.read_csv(temp_dir_path+'alw_intronic.bed',delimiter="\t",index_col=None,header=None)
        five_or_three_prime_gene_flanks_final = pd.read_csv(temp_dir_path+'five_or_three_prime_gene_flanks_final.bed',delimiter="\t",index_col=None,header=None)
        five_prime_only_gene_flanks_final = pd.read_csv(temp_dir_path+'five_prime_only_gene_flanks_final.bed',delimiter="\t",index_col=None,header=None)
        three_prime_only_gene_flanks_final = pd.read_csv(temp_dir_path+'three_prime_only_gene_flanks_final.bed',delimiter="\t",index_col=None,header=None)
        intergenic = pd.read_csv(temp_dir_path+'intergenic.bed',delimiter="\t",index_col=None,header=None)
    
        alw_exonic['segment'] = 'alw_exonic'
        alw_intronic['segment'] = 'alw_intronic'
        exonic_or_intronic['segment'] = 'exonic_or_intronic'
        intergenic['segment'] = 'intergenic'
        five_or_three_prime_gene_flanks_final['segment'] = 'five_or_three_prime_gene_flank'
        five_prime_only_gene_flanks_final['segment'] = 'five_prime_only_gene_flank'
        three_prime_only_gene_flanks_final['segment'] = 'three_prime_only_gene_flank'
        
        gene_segments = pd.concat([alw_exonic,alw_intronic,exonic_or_intronic,five_or_three_prime_gene_flanks_final,five_prime_only_gene_flanks_final,three_prime_only_gene_flanks_final,intergenic])
        gene_segments = gene_segments.sort_values([0,1,2],ascending=True).reset_index(drop=True)
        gene_segments[3] = gene_segments['segment']+';'+gene_segments[3]
        gene_segments[3] = gene_segments.groupby([0,1,2,4,5])[3].transform(lambda x: ' OR '.join(x))
        gene_segments = gene_segments[[0,1,2,3,4,5]].drop_duplicates()
    
        gene_segments[[0,1,2,3,4,5]].to_csv(temp_dir_path+'genome_segmentation_'+sel_strand+'.bed', sep=str('\t'),header=False,index=None)
    
    genome_segmentation_plus = pd.read_csv(temp_dir_path+'genome_segmentation_+.bed',delimiter="\t",index_col=None,header=None)
    genome_segmentation_minus = pd.read_csv(temp_dir_path+'genome_segmentation_-.bed',delimiter="\t",index_col=None,header=None)
    
    # Now, to get antisense regions, we overlap genic regions on one strand with intergenic on the other
    genome_segmentation_plus.loc[~genome_segmentation_plus[3].str.contains('|'.join(['intergenic','five_or_three_prime_gene_flank','five_prime_only_gene_flank','three_prime_only_gene_flank']))].to_csv(temp_dir_path+'genic_segments_plus.bed', sep=str('\t'),header=False,index=None)
    genome_segmentation_minus.loc[genome_segmentation_minus[3].str.contains('|'.join(['intergenic','five_or_three_prime_gene_flank','five_prime_only_gene_flank','three_prime_only_gene_flank']))].to_csv(temp_dir_path+'intergenic_segments_minus.bed', sep=str('\t'),header=False,index=None)
    command = 'bedtools intersect -wb -a '+temp_dir_path+'intergenic_segments_minus.bed -b '+temp_dir_path+'genic_segments_plus.bed -S > '+temp_dir_path+'antisense_segments_minus.bed'
    out = subprocess.check_output(command, shell=True)
    
    genome_segmentation_minus.loc[~genome_segmentation_minus[3].str.contains('|'.join(['intergenic','five_or_three_prime_gene_flank','five_prime_only_gene_flank','three_prime_only_gene_flank']))].to_csv(temp_dir_path+'genic_segments_minus.bed', sep=str('\t'),header=False,index=None)
    genome_segmentation_plus.loc[genome_segmentation_plus[3].str.contains('|'.join(['intergenic','five_or_three_prime_gene_flank','five_prime_only_gene_flank','three_prime_only_gene_flank']))].to_csv(temp_dir_path+'intergenic_segments_plus.bed', sep=str('\t'),header=False,index=None)
    command = 'bedtools intersect -wb -a '+temp_dir_path+'intergenic_segments_plus.bed -b '+temp_dir_path+'genic_segments_minus.bed -S > '+temp_dir_path+'antisense_segments_plus.bed'
    out = subprocess.check_output(command, shell=True)
    
    # now substract antisense regions from intergenic regions
    command = 'bedtools subtract -a '+temp_dir_path+'intergenic_segments_plus.bed -b '+temp_dir_path+'antisense_segments_plus.bed > '+temp_dir_path+'refined_intergenic_segments_plus.bed'
    out = subprocess.check_output(command, shell=True)
    command = 'bedtools subtract -a '+temp_dir_path+'intergenic_segments_minus.bed -b '+temp_dir_path+'antisense_segments_minus.bed > '+temp_dir_path+'refined_intergenic_segments_minus.bed'
    out = subprocess.check_output(command, shell=True)
    
    antisense_segments_minus = pd.read_csv(temp_dir_path+'antisense_segments_minus.bed',delimiter="\t",index_col=None,header=None)
    antisense_segments_minus[9] = 'ANTISENSE:'+antisense_segments_minus[9]
    antisense_segments_minus = antisense_segments_minus[[0,1,2,9,4,5]]
    antisense_segments_minus.columns = [0,1,2,3,4,5]
    antisense_segments_plus = pd.read_csv(temp_dir_path+'antisense_segments_plus.bed',delimiter="\t",index_col=None,header=None)
    antisense_segments_plus[9] = 'ANTISENSE:'+antisense_segments_plus[9]
    antisense_segments_plus = antisense_segments_plus[[0,1,2,9,4,5]]
    antisense_segments_plus.columns = [0,1,2,3,4,5]
    
    refined_intergenic_segments_minus = pd.read_csv(temp_dir_path+'refined_intergenic_segments_minus.bed',delimiter="\t",index_col=None,header=None)
    refined_intergenic_segments_plus = pd.read_csv(temp_dir_path+'refined_intergenic_segments_plus.bed',delimiter="\t",index_col=None,header=None)
    genic_segments_minus = genome_segmentation_minus.loc[~genome_segmentation_minus[3].str.contains('intergenic')].reset_index(drop=True)
    genic_segments_plus = genome_segmentation_plus.loc[~genome_segmentation_plus[3].str.contains('intergenic')].reset_index(drop=True)
    
    genome_segmentation = pd.concat([genic_segments_plus,refined_intergenic_segments_plus,antisense_segments_plus,
                                     genic_segments_minus,refined_intergenic_segments_minus,antisense_segments_minus]).sort_values([0,1,2],ascending=True).reset_index(drop=True)
    
    # here, marking regions as flanking has the priority over calling them antisense (for another gene)
    genome_segmentation = genome_segmentation.drop_duplicates([0,1,2,5]).sort_values([0,1,2]).reset_index(drop=True)
    genome_segmentation.to_csv(temp_dir_path+'exon_intron_gs.non_sorted.bed', sep=str('\t'),header=False,index=None)
    command = 'bedtools sort -i '+temp_dir_path+'exon_intron_gs.non_sorted.bed > '+options.output_exon_intron_gs_file_path
    out = subprocess.check_output(command, shell=True)
    
    ### NOW, get binning of genes by bin_size
    # remove intergenic and antisense segments
    sel_genome_segmentation = genome_segmentation.loc[~genome_segmentation[3].str.contains('|'.join(['ANTISENSE','intergenic']))]
    # leave only gene ids that don't have overlaps with other gene ids - NO, otherwise too few genes are left
    # bad_gene_ids = list(sel_genome_segmentation.loc[sel_genome_segmentation[3].str.contains(',')][3].str.split(';',expand=True)[1].unique())
    # bad_gene_ids = list(pd.Series((','.join(bad_gene_ids)).split(',')).unique())
    # sel_genome_segmentation = sel_genome_segmentation.loc[~sel_genome_segmentation[3].str.contains('|'.join(bad_gene_ids))]
    # # relax on alw exonic or alw intronic segments for now 
    # sel_genome_segmentation[3] = sel_genome_segmentation[3].str.split(';',expand=True)[1]
    # sel_genome_segmentation = sel_genome_segmentation.groupby([0,3,4,5]).agg({1:min,2:max}).reset_index()[[0,1,2,3,4,5]].sort_values([0,1,2,5]).reset_index(drop=True)
    
    sel_genome_segmentation.to_csv(temp_dir_path+'sel_for_binning_genome_segmentation.non_sorted.bed', sep=str('\t'),header=False,index=None)

    command = 'sort-bed '+temp_dir_path+'sel_for_binning_genome_segmentation.non_sorted.bed > '+temp_dir_path+'sel_for_binning_genome_segmentation.bedops_sorted.bed'
    out = subprocess.check_output(command, shell=True)
    
    command = 'bedops --chop '+str(bin_size)+' '+temp_dir_path+'sel_for_binning_genome_segmentation.bedops_sorted.bed | bedtools sort > '+temp_dir_path+'chopped_segments.bed'
    out = subprocess.check_output(command, shell=True)

    command = 'bedtools sort -i '+temp_dir_path+'sel_for_binning_genome_segmentation.bedops_sorted.bed > '+temp_dir_path+'sel_for_binning_genome_segmentation.bed'
    out = subprocess.check_output(command, shell=True)
    
    command = 'bedtools intersect -loj -sorted -a '+temp_dir_path+'chopped_segments.bed -b '+temp_dir_path+'sel_for_binning_genome_segmentation.bed'+\
    ' | bedtools groupby -g 1,2,3,8,9 -c 7 -o distinct | '+"""awk '{ print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}'"""+' | bedtools sort > '+options.output_binned_gs_file_path
    out = subprocess.check_output(command, shell=True)
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)