#! /usr/bin/env python
# (c) 2019. All Rights Reserved.
# Code written by: 
#       Maxim
#       Hugo
#       Sean Powell (sean.powell@imba.oeaw.ac.at)

__version__ = '0.0.1-dev'

import os
import sys
import time
import glob
import click
import tempfile

import pandas as pd
import numpy as np
import multiprocessing as mp

from snHiC.lib import coolerAnalysisShared as css
from scipy.ndimage.filters import gaussian_filter1d

COOLER_EXT = '.cool'
COOLER_FILE_NAME = '*{sample_id}*.{resolution}{cooler_ext}'
TMP_COOLER_FILE_NAME = '{prefix}.{i}.{resolution}{cooler_ext}'

MIN_LOOP_SIZE = 0
MAX_LOOP_SIZE = 1e6
SMOOTH_SIGMA = 2

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--cooler-dir', 
               type = click.Path(exists = True, resolve_path = True),
               required = True, 
               help = 'The path with the single sample cooler files.')
@click.option('--output-file',
              type = str, 
              required = True,
              help = 'The name of the output file with the bootstrap values.')
@click.option('--scaling-ref-file', 
               required = True,
               type = click.Path(exists = True, resolve_path = True),
               help = 'A merged cool file (.mcool) with the scaling that can be used as a reference.')
@click.option('--sample-ids',
              type = str, 
              required = True, 
              help = 'Comma separated list of the sample IDs.')
@click.option('--cooler-file-prefix',
              type = str, 
              required = True, 
              help = 'The prefix of the name of the temperary cooler files.')
@click.option('--tmp-dir', 
               type = click.Path(exists = True, resolve_path = True),
               required = False,
               help = 'The temporary directory for the cooler files. If none given the system tmp will be used.')
@click.option('--resolution', 
               default = 10000, 
               show_default = True,
               type = int, 
               help = 'The resolution of the cooler files.')
@click.option('--threads', 
               default = 8, 
               show_default = True,
               type = int, 
               help = 'The concurrent threads that will be started for the bootstrapping.')
@click.option('--iterations', 
               default = 1000, 
               show_default = True,
               type = int, 
               help = 'The iterations that will be bootstrapped.')
def main(cooler_dir, output_file, scaling_ref_file, sample_ids, cooler_file_prefix, tmp_dir, resolution, threads, iterations):
    '''Bootstrap the loop sizes.
    '''

    if tmp_dir == None:
        parent_tmp_dir = tempfile.gettempdir()
        tmp_dir = os.path.join(parent_tmp_dir, 'snhic')
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

    if not os.path.exists(tmp_dir):
        parent_tmp_dir = tempfile.gettempdir()
        tmp_dir = os.path.join(parent_tmp_dir, 'snhic')
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

    if threads <= 0:
        threads = 1

    sample_ids = sample_ids.split(',')
    sample_ids = [int(i) for i in sample_ids]

    results = get_results(sample_ids, resolution, scaling_ref_file, cooler_dir, cooler_file_prefix, tmp_dir, SMOOTH_SIGMA, MIN_LOOP_SIZE, MAX_LOOP_SIZE, threads, iterations)
    
    #(_, cooler_file_name) = os.path.split(cooler_file)
    #(cooler_name, _) = os.path.splitext(cooler_file_name)

    f = open(output_file, 'w')
    f.write('{}\n'.format(cooler_file_prefix))
    f.write('\n'.join([str(r[1]) for r in results]))
    f.close()

def get_results(sample_ids, resolution, scaling_ref_file, cooler_dir, cooler_prefix, tmp_dir, smooth_sigma, min_loop_size, max_loop_size, max_concurrent, jobs):
    manager = mp.Manager()
    return_list = manager.list()
    running_jobs = list()


    random_sample_ids = getSampleIDList(sample_ids, jobs)
    j = 0
    while True:
        if len(return_list) < jobs:
            for i in range(max_concurrent - len(running_jobs)):
                if j >= len(random_sample_ids):
                    break

                p = mp.Process(target = getLoopSize, args = (random_sample_ids[j], resolution, scaling_ref_file, cooler_dir, cooler_prefix, smooth_sigma, min_loop_size, max_loop_size, tmp_dir, j, return_list))
                running_jobs.append(p)
                p.start()
                j += 1

        #print('concurrent jobs: {}'.format(len(running_jobs)))

        for (i, proc) in enumerate(running_jobs):
            #proc.join()

            if not proc.is_alive():
                running_jobs.pop(i)

        #print('concurrent jobs: {}'.format(len(running_jobs)))

        if len(return_list) >= jobs:
            for running_job in running_jobs:
                running_job.terminate()

            break

        time.sleep(0.1)

    return return_list
                
def getLoopSize(sample_ids, resolution, scaling_ref_file, cooler_dir, cooler_prefix, smooth_sigma, min_loop_size, max_loop_size, tmp_dir, iteration, return_list):
    # random sample with replacement, this is basicaaly the bootstrap method
    #_sample_ids = np.random.choice(sample_ids, size = len(sample_ids), replace = True)
    
    #print(sample_ids)
    cooler_files = list()
    for sample_id in sample_ids:
        cooler_file_name = COOLER_FILE_NAME.format(sample_id = sample_id, resolution = resolution, cooler_ext = COOLER_EXT)
        cooler_file = os.path.join(cooler_dir, cooler_file_name)
       
        _cooler_files = glob.glob(cooler_file)
        if len(_cooler_files) == 0:
            return
            #return return_list

        cooler_files.append(_cooler_files[0])

    # merge the cooler files 
    merged_cooler_file_name = TMP_COOLER_FILE_NAME.format(prefix = cooler_prefix, resolution = resolution, i = iteration, cooler_ext = COOLER_EXT)
    merged_cooler_file = os.path.join(tmp_dir, merged_cooler_file_name)

    css.mergeCoolerFiles(cooler_files, merged_cooler_file)

    # get the loop size
    vals = css.getCoolerScaling(merged_cooler_file, scaling_ref_file)
    dlogx = np.diff(np.log(vals[1]))
    slope = np.diff(np.log(vals[0])) / dlogx
    mask = (vals[1][1:] * resolution > min_loop_size) & (vals[1][1:] * resolution < max_loop_size)
    mask = np.logical_not(mask)

    _slope = np.copy(slope)
    _slope = gaussian_filter1d(_slope, smooth_sigma)

    _slope[mask] = -1000
    loop_size_filt_idx = np.argmax(_slope)

    loop_size_from = vals[1][loop_size_filt_idx - 1]
    loop_size_to = vals[1][loop_size_filt_idx]

    # delete the tmp cooler file
    os.remove(merged_cooler_file)

    print('iteration {}: sample_ids: {}'.format(iteration, sample_ids))
    print('iteration {}: loop sizes: {}'.format(iteration, (loop_size_from * resolution, loop_size_to * resolution)))
    return_list.append((loop_size_from * resolution, loop_size_to * resolution))

    #return return_list 

def getSampleIDList(sample_ids, count, max_tries = None):
    if max_tries == None:
        max_tries = count * 10

    tries = 0
    sample_id_sets = dict()
    while tries < max_tries:
        tries += 1
        _sample_ids = np.random.choice(sample_ids, size = len(sample_ids), replace = True)
        _sample_ids.sort()

        _sample_ids = tuple(_sample_ids)
        if not _sample_ids in sample_id_sets:
            sample_id_sets[_sample_ids] = None

        if len(sample_id_sets.keys()) >= count:
            break

    return list(sample_id_sets.keys())


if __name__ == "__main__":
    main()

    
