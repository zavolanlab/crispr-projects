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
import csv

""" Grouping bed file with reads, weigthning multimappers and duplicated reads """

__doc__ = "Grouping bed file with reads, weigthning multimappers and duplicated reads"

parser = ArgumentParser(description=__doc__,
                        formatter_class=RawTextHelpFormatter)

parser.add_argument("--input_gtf",
                    dest="input_gtf",
                    help="path to GTF file that was used for quantification",
                    required=True,
                    metavar="FILE",)

parser.add_argument("--input_files_dir",
                    dest="input_files_dir",
                    help="path to the directory with input gct files",
                    required=True,
                    metavar="FILE")

parser.add_argument("--suffix",
                    dest="suffix",
                    help="suffix of the file names defining the quantification - gene with reads, gene with fragments, exons",
                    required=True,
                    metavar="FILE")

parser.add_argument("--temp_dir",
                dest="temp_dir",
                help="path to the directory with temporary files",
                required=True,
                metavar="FILE")

parser.add_argument("--sample_list",
                    dest="sample_list",
                    nargs="*",
                    help="names of the samples",
                    required=True)

parser.add_argument("--output_file",
                    dest="output_file",
                    help="path to output file",
                    required=True,
                    metavar="FILE",)

syserr = sys.stderr.write

def main():

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    input_gtf = options.input_gtf
    sample_list = options.sample_list
    input_files_dir = options.input_files_dir
    suffix = options.suffix
    temp_dir = options.temp_dir
    out_file_path = options.output_file
    
    gtf = pd.read_csv(input_gtf,delimiter="\t",index_col=None,header=None,skiprows=1)
    genes = gtf.loc[gtf[2]=='gene'].reset_index(drop=True)
    genes['gene_id'] = genes[8].str.split('gene_id "',expand=True)[1].str.split('";',expand=True)[0]
    genes['gene_name'] = genes[8].str.split('gene_name "',expand=True)[1].str.split('";',expand=True)[0]
    genes['gene_type'] = genes[8].str.split('gene_type "',expand=True)[1].str.split('";',expand=True)[0]
    genes = genes[[0,3,4,6,'gene_name','gene_type','gene_id']].drop_duplicates('gene_id').reset_index(drop=True)
    genes.index = genes['gene_id']
    genes.columns = ['chr','start','end','strand','gene_name','gene_type','gene_id']
    
    genes[['chr','strand','gene_type']] = genes[['chr','strand','gene_type']].astype('category')
    
    # find all relevant files and match them with samples
    command = """mkdir -p """+temp_dir+"""; find """+input_files_dir+""" -name '*"""+suffix+"""' > """+temp_dir+"""/rnaseqc_file_paths.tsv"""
    out = subprocess.check_output(command, shell=True)
    rnaseqc_file_paths = pd.read_csv(temp_dir+"""/rnaseqc_file_paths.tsv""",delimiter="\t",index_col=None,header=None)
    rnaseqc_file_paths = rnaseqc_file_paths.loc[rnaseqc_file_paths[0].str.contains('|'.join(sample_list))].reset_index(drop=True)
    a = []
    for elem in list(rnaseqc_file_paths[0]):
        found=False
        for sample in sample_list:
            if sample in elem:
                a.append(sample)
                found = True
                break
        if not found:
            a.append("")
    rnaseqc_file_paths['sample'] = a
    rnaseqc_file_paths = rnaseqc_file_paths.loc[rnaseqc_file_paths['sample']!=""].reset_index(drop=True)
    rnaseqc_file_paths = rnaseqc_file_paths.rename(columns = {0:'path'})
    
    order_df = genes[['gene_id']].copy().reset_index(drop=True)
    num_of_samples = len(rnaseqc_file_paths)
    i = 1
    a = []
    b = []
    for index, row in rnaseqc_file_paths.iterrows():
        tmp = pd.read_csv(row['path'],delimiter="\t",index_col=None,header=0,usecols = [0,2])
        tmp.columns = ['gene_id','r']
        a.append(pd.merge(order_df,tmp,how='left',on=['gene_id'])['r'])
        b.append(row['sample'])
        if i%50==0:
            print('[INFO] '+str(i)+' out of '+str(num_of_samples)+' loaded')
        i=i+1
    counts = pd.DataFrame(np.transpose(a),columns = b)
    genes = pd.concat([genes.reset_index(drop=True),counts],axis=1)
    genes[list(genes.columns)[7:]] = genes[list(genes.columns)[7:]].fillna(0).astype('int')
    
    genes.to_csv(out_file_path, sep=str('\t'),header=True,index=None,quoting=csv.QUOTE_NONE)
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)