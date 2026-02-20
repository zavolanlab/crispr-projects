#!/usr/bin/env python3

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess
import pandas as pd
import numbers
import time

def main():
    """ 
    """

    __doc__ = ""

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument("--max_read_lengths_aggr",
                        dest="max_read_lengths_aggr_path",
                        help="Path to the max_read_lengths_aggr file",
                        required=True,
                        metavar="FILE",)
    
    parser.add_argument("--cluster_log_dir",
                        dest="cluster_log_dir",
                        help="directory of cluster logs",
                        required=True,
                        metavar="DIR")
    
    parser.add_argument("--threads_number",
                        dest="threads_number",
                        help="number of threads per run",
                        required=True,
                        metavar="INT")    
    
    parser.add_argument("--docker_image",
                        dest="docker_image",
                        help="URL for docker image to pull",
                        required=True,
                        metavar="URL")
    
    parser.add_argument("--docker_image_target_path",
                        dest="docker_image_target_path",
                        help="path for docker image to download to",
                        required=True,
                        metavar="FILE")    

    parser.add_argument("--results_dir",
                        dest="results_dir",
                        help="directory where indices will be stored",
                        required=True,
                        metavar="DIR")    
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    # create the results directory and scripts subdirectory within it
    out = subprocess.check_output('mkdir -p '+options.results_dir, shell=True)
    out = subprocess.check_output('mkdir -p '+os.path.join(options.results_dir,'scripts'), shell=True)
    
    # pull docker image
    out = subprocess.check_output('singularity pull -F '+options.docker_image_target_path+' '+options.docker_image, shell=True)
    
    tasks_df = pd.read_csv(options.max_read_lengths_aggr_path,delimiter="\t",index_col=None,header=0)
    tasks_df = tasks_df.drop_duplicates(['organism','max_read_length_estimate']).reset_index(drop=True) # ensure no duplicated jobs are submitted
    # loop over tasks and submit slurm jobs
    JobIDs = []
    for index,row in tasks_df.iterrows():
        cur_indices_dir = os.path.join(options.results_dir,row['organism']+'_'+str(row['max_read_length_estimate']))
        # create subdirectory for this task
        out = subprocess.check_output('mkdir -p '+cur_indices_dir, shell=True)
        # create slurm script
        slurm_script_path = os.path.join(options.results_dir,'scripts',row['organism']+'_'+str(row['max_read_length_estimate'])+'.sbatch')
        
        dirs_to_bind = ','.join(list(pd.Series([os.path.dirname(row['genome_file']),os.path.dirname(row['gtf_file']),cur_indices_dir]).unique()))
        f = open(slurm_script_path, "w")
        preambula = \
"""#!/bin/bash

#SBATCH --job-name=create_star_index_"""+row['organism']+'_'+str(row['max_read_length_estimate'])+"""
#SBATCH --time=5:00:00
#SBATCH --qos=6hours
#SBATCH --output="""+os.path.join(options.cluster_log_dir,'create_STAR_indices.'+row['organism']+'_'+str(row['max_read_length_estimate'])+'.out')+"""
#SBATCH --error="""+os.path.join(options.cluster_log_dir,'create_STAR_indices.'+row['organism']+'_'+str(row['max_read_length_estimate'])+'.err')+"""
#SBATCH --cpus-per-task="""+str(options.threads_number)+"""
#SBATCH --mem=50G

source ~/.bashrc"""
        f.write(preambula+'\n')                         
        command = \
"""
singularity exec --bind """+dirs_to_bind+' '+options.docker_image_target_path+""" \
STAR --runThreadN """+str(options.threads_number)+""" \
--runMode genomeGenerate \
--genomeDir """+cur_indices_dir+""" \
--outFileNamePrefix """+os.path.join(cur_indices_dir,'STAR_')+""" \
--genomeFastaFiles """+row['genome_file']+""" \
--sjdbGTFfile """+row['gtf_file']+""" \
--sjdbOverhang """+str(row['max_read_length_estimate']-1)
        f.write(command)
        f.close()
        out = subprocess.check_output('sbatch '+slurm_script_path, shell=True)
        JobID = out.decode().replace('Submitted batch job ','').replace('\n','')
        JobIDs.append(JobID)
    
    # check status every minute, and do not run the job longer than 12 hours
    finished = False
    statuses = ['NOT_COMPLETED']
    start_time = time.time()
    while ((time.time()-start_time)<60*60*12):
        finished, statuses = check_jobs_status(JobIDs)
        if finished:
            break
        else:
            time.sleep(60)
    
    # if all jobs has status COMPLETED, than proceed, else report error
    report = [status for status in statuses if status!='COMPLETED']
    if len(report)==0:
        # create a flag that all tasks were successfully executed
        out = subprocess.check_output('echo -n "success" > '+os.path.join(options.results_dir,'success.txt'), shell=True)    
    else:
        sys.stderr.write("Not all jobs were successfully completed")
        sys.exit(1)
        
def check_jobs_status(JobIDs):
    finished = True
    statuses = []
    for JobIDs in JobIDs:
        out = subprocess.check_output('sacct --job '+str(JobIDs)+' -o JobID,JobName,State ', shell=True)
        status = out.decode().split('\n')[2].strip().split(' ')[-1]
        if not (status in ['COMPLETED','CANCELLED','BOOT_FAIL','DEADLINE','FAILED','NODE_FAIL','OUT_OF_MEMORY','PREEMPTED','TIMEOUT']):
            finished = False
        statuses.append(status)
    return finished, statuses    
    
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)