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
import itertools

def main():
    """ get cigar stats from sam file"""

    __doc__ = " get cigar stats from sam file"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--category_specific_sam_file",
                        dest="category_specific_sam_file",
                        help="Path to the input category_specific_sam_file",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--category_specific_read_names",
                        dest="category_specific_read_names",
                        help="path for the input category_specific_read_names",
                        required=True,
                        metavar="FILE")
    parser.add_argument("--output_file",
                        dest="output_file",
                        help="path for output file cigar_stats.tsv",
                        required=True,
                        metavar="FILE")
  
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    category_specific_sam_file = options.category_specific_sam_file
    category_specific_read_names = options.category_specific_read_names
    output_file = options.output_file
    
    category_specific_sam_file = pd.read_csv(category_specific_sam_file,delimiter="\t",index_col=None,header=None,usecols=[0,5,11]) # only read name and weight
    category_specific_sam_file.columns = ['read_name','cigar','MD']
    category_specific_sam_file['MD'] = category_specific_sam_file['MD'].str.replace('MD:Z:','')
    category_specific_sam_file['t']=1
    category_specific_sam_file = pd.merge(category_specific_sam_file.drop('t',1),category_specific_sam_file.groupby('read_name').agg({'t':sum}).reset_index(),how='inner',on=['read_name'])
    category_specific_sam_file['w'] = 1/category_specific_sam_file['t']
    gr = category_specific_sam_file.groupby(['cigar','MD']).agg({'w':sum}).reset_index().sort_values('w',ascending=False).reset_index(drop=True)
    gr.to_csv(output_file, sep=str('\t'),header=True,index=None)
            
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)