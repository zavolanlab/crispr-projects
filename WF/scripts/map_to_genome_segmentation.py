#!/usr/bin/env python3

# -----------------------------------------------------------------------------
# Author : Aleksei Mironov
# Company: Mihaela Zavolan, Biozentrum, Basel
# This script is part of the Zavolan lab ZARP pipeline.
# -----------------------------------------------------------------------------

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd

def main():
    """Mapping to the genome segmentation and getting counts and coverage"""

    __doc__ = "Mapping to the genome segmentation and getting counts and coverage"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--grouped_bed",
                        dest="grouped_bed_file_path",
                        help="Path to the input grouped bed file produced with groupBamToBedFile",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--genome_segmentation",
                        dest="genome_segmentation_file_path",
                        help="path for the genome segmentation bed file",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--wd_counts",
                        dest="wd_counts_file_path",
                        help="path for the output file with counts per genomic segment",
                        required=True,
                        metavar="FILE")    
    parser.add_argument("--wd_coverage",
                        dest="wd_coverage_file_path",
                        help="path for the output file with coverage histograms per genomic segment",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--directionality",
                        dest="directionality",
                        help="Library directionality (first read - reverse (SR, -S) or first read - forward (SF, -s))",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--temp_dir",
                        dest="temp_dir",
                        help="path where to put temporary files",
                        required=True,
                        metavar="FILE")
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    genome_segmentation_file = options.genome_segmentation_file_path
    grouped_bed_file_path = options.grouped_bed_file_path
    wd_counts_output_file = options.wd_counts_file_path
    wd_coverage_output_file = options.wd_coverage_file_path
    s_param = '-'+options.directionality
    prefix = options.temp_dir+'/'
    
    output_dir = '/'.join((wd_counts_output_file.split('/')[:-1]))+'/'

    # create output directory if they don't exist
    out = subprocess.check_output('mkdir -p '+output_dir, shell=True)
    out = subprocess.check_output('mkdir -p '+'/'.join((wd_coverage_output_file.split('/')[:-1]))+'/', shell=True)

    # check that bed file is not empty, otherwise create empty output files
    if not (os.path.isfile(options.grouped_bed_file_path) and (os.stat(options.grouped_bed_file_path).st_size > 0)):
        out = subprocess.check_output('touch '+wd_counts_output_file, shell=True)
        out = subprocess.check_output('touch '+wd_coverage_output_file, shell=True)
        
    bed_file = pd.read_csv(options.grouped_bed_file_path,delimiter="\t",index_col=None,header=None)
    bed_file.columns = [0,1,2,'name','w',5,'d','wd','d_cat']
    
    d_cats = list(bed_file['d_cat'].unique())
    d_cats.sort()

    a = [] # list of file names
    to_rm = [] # list of temporary files
    for mm_mode in ['uniquely_mapped','multimapped']:
        if mm_mode=='uniquely_mapped':
            cur_bed_file = bed_file.loc[bed_file['w']==1].reset_index(drop=True)
        else:
            cur_bed_file = bed_file.loc[bed_file['w']<1].reset_index(drop=True)
        
        for dupl_category in d_cats:
            d_cur_bed_file = cur_bed_file.loc[cur_bed_file['d_cat']==dupl_category].reset_index(drop=True)
            if len(d_cur_bed_file)>0:
                non_sorted_file_name_in = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.nonsorted.bed'
                file_name_in = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.bed'
                d_cur_bed_file[[0,1,2,'name','wd',5]].to_csv(non_sorted_file_name_in, sep=str('\t'),header=False,index=None) # with counts of duplicates
                command = 'bedtools sort -i '+non_sorted_file_name_in+' > '+file_name_in
                out = subprocess.check_output(command, shell=True)
                to_rm.append(non_sorted_file_name_in)
                to_rm.append(file_name_in)
                file_name_out = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.intersected.bed'
                # account for reads that span multiple features
                out = subprocess.check_output('bedtools intersect -sorted -a '+file_name_in+' -b '+genome_segmentation_file+' -c '+s_param+' > '+file_name_out, shell=True)
                if out.decode()=='': # if successful
                    to_rm.append(file_name_out)
                    d_cur_bed_file = pd.read_csv(file_name_out,delimiter="\t",index_col=None,header=None)
                    d_cur_bed_file = d_cur_bed_file.loc[d_cur_bed_file[6]>0].reset_index(drop=True)
                    if len(d_cur_bed_file)>0:
                        d_cur_bed_file['MF_W'] = 1/d_cur_bed_file[6] # weight for spanning multiple features
                        d_cur_bed_file[4] = (d_cur_bed_file[4]*d_cur_bed_file['MF_W']*1000).astype('int') # scale to 1000
                        d_cur_bed_file[list(range(0,6))].to_csv(file_name_out, sep=str('\t'),header=False,index=None) # rewrite the output with weights for spanning multiple features

                        file_name_out_positive = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.intersected.+.bed'
                        d_cur_bed_file.loc[d_cur_bed_file[5]=='+'][list(range(0,6))].to_csv(file_name_out_positive, sep=str('\t'),header=False,index=None)
                        to_rm.append(file_name_out_positive)
                        file_name_out_negative= prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.intersected.-.bed'
                        d_cur_bed_file.loc[d_cur_bed_file[5]=='-'][list(range(0,6))].to_csv(file_name_out_negative, sep=str('\t'),header=False,index=None)
                        to_rm.append(file_name_out_negative)
                        file_name_in = file_name_out
                        ### wd sum
                        wd_file_name_out = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.wd_sum.bed'
                        out = subprocess.check_output('bedtools map -a '+genome_segmentation_file+' -b '+file_name_in+' -null 0 '+s_param+' -c 5 -o sum > '+wd_file_name_out, shell=True)
                        to_rm.append(wd_file_name_out)
                        ### get coverage summary

                        file_name_out_positive_partitions = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.partitions.+.bed'
                        file_name_out_negative_partitions = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.partitions.-.bed'

                        file_name_out_positive_coverage = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.coverage.+.bed'
                        file_name_out_negative_coverage = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.coverage.-.bed'
                        file_name_out_coverage = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.coverage.bed'

                        coverage_intersect_file_name_out = prefix+'.wd.'+mm_mode+'.'+str(dupl_category)+'.coverage_intersect.bed'
                        if len(d_cur_bed_file.loc[d_cur_bed_file[5]=='+'])>0:

                            # bedops requires files to be sorted with sort-bed
                            command = 'sort-bed '+file_name_out_positive+' > '+file_name_out_positive+'.sorted && rm '+file_name_out_positive+' && mv '+file_name_out_positive+'.sorted '+file_name_out_positive
                            out = subprocess.check_output(command, shell=True)
                            
                            command = 'bedops --partition '+file_name_out_positive+' > '+file_name_out_positive_partitions
                            out = subprocess.check_output(command, shell=True)
                            
                            to_rm.append(file_name_out_positive_partitions)
                            command = 'bedtools map -a '+file_name_out_positive_partitions+' -b '+file_name_out_positive+' -null 0 -c 5 -o sum > '+file_name_out_positive_coverage
                            out = subprocess.check_output(command, shell=True)
                            to_rm.append(file_name_out_positive_coverage)
                            command = 'sed -i "s/$/\t0\t+/" '+file_name_out_positive_coverage
                            out = subprocess.check_output(command, shell=True)
                        else:
                            command = 'touch '+file_name_out_positive_coverage # create empty file
                            out = subprocess.check_output(command, shell=True)
                            to_rm.append(file_name_out_positive_coverage)

                        if len(d_cur_bed_file.loc[d_cur_bed_file[5]=='-'])>0:

                            # bedops requires files to be sorted with sort-bed
                            command = 'sort-bed '+file_name_out_negative+' > '+file_name_out_negative+'.sorted && rm '+file_name_out_negative+' && mv '+file_name_out_negative+'.sorted '+file_name_out_negative
                            out = subprocess.check_output(command, shell=True)                            
                            
                            command = 'bedops --partition '+file_name_out_negative+' > '+file_name_out_negative_partitions
                            out = subprocess.check_output(command, shell=True)
                            
                            to_rm.append(file_name_out_negative_partitions)
                            command = 'bedtools map -a '+file_name_out_negative_partitions+' -b '+file_name_out_negative+' -null 0 -c 5 -o sum > '+file_name_out_negative_coverage
                            out = subprocess.check_output(command, shell=True)
                            to_rm.append(file_name_out_negative_coverage)

                            command = 'sed -i "s/$/\t0\t-/" '+file_name_out_negative_coverage
                            out = subprocess.check_output(command, shell=True)
                        else:
                            command = 'touch '+file_name_out_negative_coverage # create empty file
                            out = subprocess.check_output(command, shell=True)
                            to_rm.append(file_name_out_negative_coverage)

                        command = 'cat '+file_name_out_negative_coverage+' '+file_name_out_positive_coverage+' | bedtools sort > '+file_name_out_coverage
                        out = subprocess.check_output(command, shell=True)
                        to_rm.append(file_name_out_coverage)
                        command = 'bedtools intersect -sorted -a '+genome_segmentation_file+' -b '+file_name_out_coverage+' -wo '+s_param+' | cut -f1-4,6,10,13 > '+coverage_intersect_file_name_out
                        out = subprocess.check_output(command, shell=True)
                        to_rm.append(coverage_intersect_file_name_out)
                        coverage_intersect = pd.read_csv(coverage_intersect_file_name_out,delimiter="\t",index_col=None,header=None)
                        coverage_intersect = coverage_intersect.rename(columns={5:'r',6:'overlap'}).groupby([0,1,2,3,4,'r']).agg({'overlap':sum}).reset_index()
                        coverage_intersect = coverage_intersect.sort_values([0,1,2,3,'r']).reset_index(drop=True)
                        coverage_intersect['score'] = dupl_category
                        coverage_intersect['r'] = coverage_intersect['r']/1000
                        coverage_intersect['mm_mode'] = ('um' if mm_mode=='uniquely_mapped' else 'mm')
                        coverage_intersect[[0,1,2,3,'score',4,'r','overlap','mm_mode']].to_csv(coverage_intersect_file_name_out, sep=str('\t'),header=False,index=None)
                        a.append([mm_mode+';'+str(dupl_category),wd_file_name_out,coverage_intersect_file_name_out])
                    else:
                        a.append([mm_mode+';'+str(dupl_category),None,None])
            else:
                a.append([mm_mode+';'+str(dupl_category),None,None])

    # merging wd counts and coverage
    res = pd.read_csv(genome_segmentation_file,delimiter="\t",index_col=None,header=None)
    res = res[[0,1,2,5,3]]
    concat_string = ''
    for elem in a:
        # wd counts
        if elem[1]!=None:
            tmp = pd.read_csv(elem[1],delimiter="\t",index_col=None,header=None)
            tmp[6] = tmp[6]/1000
            res = pd.merge(res,tmp[[0,1,2,5,3,6]].rename(columns={6:elem[0]}),how='left',on=[0,1,2,5,3])
        else:
            res[elem[0]]=0
        if elem[2]!=None:
            concat_string = concat_string+' '+elem[2]
    
    cols = list(res.columns)
    cols = ['chr','start','end','strand','segment_name']+cols[5:]
    res.columns = cols
    # res = res.loc[res[cols[5:]].sum(1)>0].reset_index(drop=True) # remove the lines with all zeros
    res.drop_duplicates(['chr','start','end',
                         'strand','segment_name']).sort_values(['chr','start','end',
                                                                'strand','segment_name']).to_csv(wd_counts_output_file, sep=str('\t'),header=True,index=None)

    # concatenate all coverage files into the final output file
    command = 'cat'+concat_string+' > '+wd_coverage_output_file
    out = subprocess.check_output(command, shell=True)
    # now we can delete all the created unnecessary files
    to_rm = [f for f in to_rm if not(f in [wd_counts_output_file, wd_coverage_output_file])]
    out = subprocess.check_output('rm '+' '.join(to_rm), shell=True)
        
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)