#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd
import numbers

def main():
    """ Adding max read length estimates and generating all unique combinations organism x read length
        for the generation of STAR indices 
    """

    __doc__ = "Adding max read length estimates and generating all unique combinations organism x read length"

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_samples",
                        dest="input_samples_file_path",
                        help="Path to the input samples file",
                        required=True,
                        metavar="FILE",)

    parser.add_argument("--results_dir",
                        dest="results_dir",
                        help="path for the results directory in which subdirectory <samples> should exist",
                        required=True,
                        metavar="DIR")

    parser.add_argument("--out_file",
                        dest="out_file_path",
                        help="path for the output file",
                        required=True,
                        metavar="FILE")    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    input_samples = pd.read_csv(options.input_samples_file_path,delimiter="\t",index_col=0,header=0)
    # expected path for the read length estimate
    input_samples['exp_read_length_estimate_file_path'] = options.results_dir+"samples/"+input_samples.index+"/read_length/"+input_samples.index+".max_read_length.txt"
    
    a = []
    for index, row in input_samples.iterrows():
        if (os.path.isfile(row['exp_read_length_estimate_file_path']) and (os.stat(row['exp_read_length_estimate_file_path']).st_size > 0)):
            exp_read_length_estimate = pd.read_csv(row['exp_read_length_estimate_file_path'],delimiter="\t",index_col=None,header=None)
            estimate = exp_read_length_estimate.iloc[0][0]
            if isinstance(estimate, numbers.Number):
                a.append(estimate)
            else:
                sys.stderr.write("file "+row['exp_read_length_estimate_file_path']+" contains a non-integer value. Exiting.")
                sys.exit(1)
        else:
            sys.stderr.write("file "+row['exp_read_length_estimate_file_path']+" does not exist or empty. Exiting.")
            sys.exit(1)
    
    input_samples['max_read_length_estimate'] = a
    res = input_samples[['organism','max_read_length_estimate','genome_file','gtf_file']].drop_duplicates().reset_index(drop=True)
    res.to_csv(options.out_file_path, sep=str('\t'),header=True,index=None)
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)