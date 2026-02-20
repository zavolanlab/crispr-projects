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
import numpy as np

def main():
    """ Grouping bed file with reads, weigthning multimappers and duplicated reads """

    __doc__ = "Grouping bed file with reads, weigthning multimappers and duplicated reads"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_bed",
                        dest="input_bed_file_path",
                        help="Path to the input bed file produced with BamToBed",
                        required=True,
                        metavar="FILE",)

    parser.add_argument("--grouped_bed",
                        dest="grouped_bed_file_path",
                        help="path for the grouped output file",
                        required=True,
                        metavar="FILE")

    parser.add_argument("--duplicates_summary",
                        dest="duplicates_summary_file_path",
                        help="path for the duplicates summary output file",
                        required=True,
                        metavar="FILE")    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    bed_file = pd.read_csv(options.input_bed_file_path,delimiter="\t",index_col=None,header=None)
    bed_file['t']=1
    gr = bed_file.groupby(3).agg({'t':sum}).reset_index()
    bed_file = pd.merge(bed_file.drop('t',1),gr,how='left',on=3)
    bed_file['w'] = 1/bed_file['t'] # weight multi-mapppers uniformely
    bed_file['d']=1

    bed_file = bed_file.groupby([0,1,2,5]).agg({'w':np.median,'d':sum}).reset_index()
    bed_file['name'] = bed_file.index
    bed_file['wd'] = bed_file['w']*bed_file['d']
    # here we define duplication categories # TO DO: allow parametrization
    bed_file['d_cat'] = pd.cut(bed_file['d'],bins=[0,10,100,1000,10000,100000,max(1000000,bed_file['d'].max())],labels=False)
    
    # create output directory if not exists
    outdir = '/'.join((options.grouped_bed_file_path.split('/')[:-1]))+'/'
    out = subprocess.check_output('mkdir -p '+outdir, shell=True)

    output_dir = os.path.dirname(options.grouped_bed_file_path)
    temporary_file = output_dir+'/'+os.path.basename(options.grouped_bed_file_path).replace('.bed','.nonsorted.bed')
    
    bed_file[[0,1,2,'name','w',5,'d','wd','d_cat']].to_csv(temporary_file, sep=str('\t'),header=False,index=None)
    
    command = 'bedtools sort -i '+temporary_file+' > '+options.grouped_bed_file_path+' && rm '+temporary_file
    out = subprocess.check_output(command, shell=True)
    
    bed_file['t'] = 1
    gr = bed_file.groupby('d_cat').agg({'t':sum,'w':np.median}).reset_index()
    gr = pd.merge(gr,bed_file.loc[bed_file['w']==1].groupby('d_cat').agg({'wd':sum}).reset_index().rename(columns={'wd':'wd_um'}),how='left',on='d_cat')
    gr = pd.merge(gr,bed_file.loc[bed_file['w']<1].groupby('d_cat').agg({'wd':sum}).reset_index().rename(columns={'wd':'wd_mm'}),how='left',on='d_cat')
    gr['prop_um'] = np.round(gr['wd_um']/gr['wd_um'].sum()*100,1) # proportion of uniquely-mapped reads in a duplication category
    gr['prop_mm'] = np.round(gr['wd_mm']/gr['wd_mm'].sum()*100,1) # proportion of multi-mapped reads in a duplication category
    gr['w'] = (1/gr['w']).astype('int') # median number of reads per loci in this duplication category (1 - unique, >1 - multi-mappers)
    
    # create output directory if not exists
    outdir = '/'.join((options.duplicates_summary_file_path.split('/')[:-1]))+'/'
    out = subprocess.check_output('mkdir -p '+outdir, shell=True)
    
    gr.to_csv(options.duplicates_summary_file_path, sep=str('\t'),header=True,index=None) # with counts of duplicates    
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)