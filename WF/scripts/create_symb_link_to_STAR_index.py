#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd
import numbers

def main():
    """  
    """

    __doc__ = ""

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--input_samples",
                        dest="input_samples_file_path",
                        help="Path to the input samples file",
                        required=True,
                        metavar="FILE",)

    parser.add_argument("--sample_id",
                        dest="sample_id",
                        help="sample_id",
                        required=True,
                        metavar="STR",)

    parser.add_argument("--max_read_length_estimate",
                        dest="max_read_length_estimate_file_path",
                        help="Path to the file where max_read_length_estimate is stored",
                        required=True,
                        metavar="FILE",)    

    parser.add_argument("--STAR_indices_dir",
                        dest="STAR_indices_dir",
                        help="directory in which STAR indices were created",
                        required=True,
                        metavar="DIR")

    parser.add_argument("--link_to_create",
                        dest="link_to_create",
                        help="symbolic link to create",
                        required=True,
                        metavar="FILE")    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    # loading input and checking consistency
    input_samples = pd.read_csv(options.input_samples_file_path,delimiter="\t",index_col=0,header=0)
    
    exp_read_length_estimate = pd.read_csv(options.max_read_length_estimate_file_path,delimiter="\t",index_col=None,header=None)
    estimate = exp_read_length_estimate.iloc[0][0]
    if not isinstance(estimate, numbers.Number):
        sys.stderr.write("file "+options.max_read_length_estimate_file_path+" contains a non-integer value. Exiting.")
        sys.exit(1)
        
    organism = input_samples.loc[options.sample_id]['organism']
    
    # expected directory where STAR index should have been created
    expected_STAR_index_dir = os.path.join(options.STAR_indices_dir,organism+'_'+str(estimate))
    
    # check that the expected dir exists and there is expected content there
    if os.path.isdir(expected_STAR_index_dir) and os.path.isfile(os.path.join(expected_STAR_index_dir,'SAindex')) and (os.stat(os.path.join(expected_STAR_index_dir,'SAindex')).st_size > 10**9):
        command = 'ln -s '+expected_STAR_index_dir+' '+options.link_to_create
        out = subprocess.check_output(command, shell=True)
    else:
        sys.stderr.write("directory "+expected_STAR_index_dir+" does not exist or contains SAindex of unexpectedly small size. Exiting.")
        sys.exit(1)
        
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)