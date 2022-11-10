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

from scipy import stats

from pandas import ExcelWriter
from mirnylib.genome import Genome
from snHiC.lib import coolerAnalysisShared as css

COOLER_EXT = '.cool'
COOLER_FILE_NAME = '*{sample_id}*.{resolution}{cooler_ext}'
TMP_COOLER_FILE_NAME = '{prefix}.{i}.{resolution}{cooler_ext}'
MAX_EXP_LOOP_SIZE = 1e6
SMOOTH_SIGMA = 2

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--condition-a', 
               type = str,
               required = True)
@click.option('--condition-b', 
               type = str,
               required = True)
@click.option('--sample-size-a', 
               type = int,
               required = True,
               help = 'The name of the first sample condition.')
@click.option('--sample-size-b', 
               type = int,
               required = True,
               help = 'The name of the second sample condition.')
@click.option('--sample-ids', 
               type = str,
               required = True,
               help = 'A comma delimited list of the sample IDs.')
@click.option('--output-file',
              type = str, 
              required = True,
              help = 'The name of the outfile in .xlsx format.')
@click.option('--genome-path', 
               required = True, 
               type = click.Path(exists = True, resolve_path = True), 
               help = 'The directory with the chomosome sequences in fasta format.')
@click.option('--ref-cooler', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The reference cooler used to filter out noise.')
@click.option('--chroms-file', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'A file with the chromosome sizes')
@click.option('--smooth-sigma', 
               default = SMOOTH_SIGMA,
               type = str, 
               show_default = True,
               help = 'The smooth sigma value for the gaussian filter.')
@click.option('--max-loop-size', 
               default = MAX_EXP_LOOP_SIZE,
               type = str,
               show_default = True,
               help = 'The maximum loop size taken into account.')
@click.option('--cooler-dir', 
               type = click.Path(exists = True, resolve_path = True),
               required = True,
               help = 'The  directory with the cooler files.')
@click.option('--cooler-ext', 
               default = COOLER_EXT,
               show_default = True,
               type = str, 
               help = 'The cooler file extension.')
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
               help = 'The concurrent threads that will be started for the random sampling.')
@click.option('--iterations', 
               default = 100, 
               show_default = True,
               type = int, 
               help = 'The iterations of the rampldom sampling.')
def main(condition_a, condition_b, sample_size_a, sample_size_b, sample_ids, output_file, genome_path, loops_file, cooler_dir, cooler_ext, tmp_dir, resolution, threads, iterations):
    '''Randomize the aggregate loops analysis.
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

    gen = Genome(genome_path, readChrms = ["#", "X"])
    loops_df = pd.read_csv(loops_file, sep = "\t")

    if threads <= 0:
        threads = 1

    if threads > iterations:
        threads = iterations

    print(condition_a)
    print(condition_b)
    sample_ids = sample_ids.split(',')
    sample_ids = [int(i) for i in sample_ids]
    sample_ids = np.array(sample_ids)

    results = get_results(loops_df, gen, condition_a, condition_b, sample_size_a, sample_size_b, sample_ids, cooler_dir, cooler_ext, tmp_dir, resolution, threads, iterations)
    #print(results)
    results_a = list()
    results_b = list()

    for (result_a, result_b) in results:
        results_a.append(result_a)
        results_b.append(result_b)
    #
    #if os.path.exists(output_file):
    # return 

    column_names = [condition_a, condition_b]

    data = np.array([np.array(results_a), np.array(results_b)])

    print(stats.wilcoxon(np.array(results_a), np.array(results_b)))

    #I = pd.Index(chr_names, name = 'iteration')
    C = pd.Index(column_names)
    df = pd.DataFrame(data.transpose(), columns = C)

    writer = ExcelWriter(output_file)
    df.to_excel(writer)
    writer.save()

def get_results(loops_df, gen, condition_a, condition_b, sample_size_a, sample_size_b, sample_ids, cooler_dir, cooler_ext, tmp_dir, resolution, max_concurrent, jobs):
    manager = mp.Manager()
    return_list = manager.list()
    running_jobs = list()

    j = 0
    while True:
        if len(return_list) < jobs:
            for i in range(max_concurrent - len(running_jobs)):
                j += 1
                if j > jobs:
                    break

                p = mp.Process(target = getAggrLoopStrength, args = (loops_df, gen, condition_a, condition_b, sample_size_a, sample_size_b, sample_ids, cooler_dir, cooler_ext, tmp_dir, resolution, j, return_list))
                running_jobs.append(p)
                p.start()

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

def getAggrLoopStrength(loops_df, gen, condition_a, condition_b, sample_size_a, sample_size_b, sample_ids, cooler_dir, cooler_ext, tmp_dir, resolution, j, return_list):
    # random sample with replacement, this is basicaaly the bootstrap method
    #loops_df = loops_df.sample(len(loops_df), replace = True, random_state = np.random.RandomState())
    #print(loops_df)
    (sample_ids_a, sample_ids_b) = randomized_and_split(sample_ids, sample_size_a, sample_size_b)
    print('iteration: {}'.format(j))
    #print('samples_a: {}'.format(sample_ids_a))
    #print('samples_b: {}'.format(sample_ids_b))
    cooler_file_a = create_tmp_cooler_file(condition_a, j, tmp_dir, sample_ids, resolution, cooler_dir, cooler_ext)
    cooler_file_b = create_tmp_cooler_file(condition_b, j, tmp_dir, sample_ids, resolution, cooler_dir, cooler_ext)

    ################
    ##### Continue here

    #print('getLoops')
    dictLoops = css.getLoops(gen, loops_df)
    #print('doPileUp')
    stored_maps_a = doPileUp(cooler_file_a, dictLoops)
    stored_maps_b = doPileUp(cooler_file_b, dictLoops)
    #print('showLoopPileUp')
    (M_a, C_a, m_count_a, c_count_a) = css.showLoopPileUp(stored_maps_a)
    (M_b, C_b, m_count_b, c_count_b) = css.showLoopPileUp(stored_maps_b)
    #print('getLoopStrength')
    strength_a = getLoopStrength(M_a, m_count_a, C_a, c_count_a)
    strength_b = getLoopStrength(M_b, m_count_b, C_b, c_count_b)

    #remove the temp cooler files
    os.remove(cooler_file_a)
    os.remove(cooler_file_b)
    
    return_list.append((strength_a, strength_b))

def randomized_and_split(sample_ids, sample_size_a, sample_size_b):
    if not len(sample_ids) == sample_size_a + sample_size_b:
        raise Exception

    np.random.seed()
    sample_ids_a = np.random.choice(sample_ids, size = sample_size_a, replace = False)
    sample_ids_b = list()
    for sample_id in sample_ids:
        if not sample_id in sample_ids_a:
            sample_ids_b.append(sample_id)

    return (sample_ids_a, sample_ids_b)

def create_tmp_cooler_file(condition_name, iteration, tmp_dir, sample_ids, resolution, cooler_dir, cooler_ext):
    cooler_files = list()
    for sample_id in sample_ids:
        cooler_file_name = COOLER_FILE_NAME.format(sample_id = sample_id, resolution = resolution, cooler_ext = cooler_ext)
        cooler_file = os.path.join(cooler_dir, cooler_file_name)
       
        _cooler_files = glob.glob(cooler_file)
        if len(_cooler_files) == 0:
            return
            #return return_list

        cooler_files.append(_cooler_files[0])

    # merge the cooler files 
    merged_cooler_file_name = TMP_COOLER_FILE_NAME.format(prefix = condition_name, resolution = resolution, i = iteration, cooler_ext = COOLER_EXT)
    merged_cooler_file = os.path.join(tmp_dir, merged_cooler_file_name)

    css.mergeCoolerFiles(cooler_files, merged_cooler_file)

    return merged_cooler_file

def getLoopStrength(M, m_count, C, c_count):
    if m_count == 0 or c_count == 0:
        return 0
    
    return css.getLoopStrength(M / m_count, C / c_count)

def doPileUp(cooler_file, dictLoops):
    (_, file_name) = os.path.split(cooler_file)
    
    separations = [100000, 150000, 250000, 500000]
    #padsize = 10 # 120 kb if resolution = 10000
    dictLoops['sep'] = dictLoops['ends'] - dictLoops['starts'] 
    

    stored_maps = dict()
    for (_, lo, hi) in zip(range(len(separations[0:-1])), separations[0:-1], separations[1:]):

        goodLoops = (dictLoops['sep'] >= lo) & (dictLoops['sep'] < hi)
        these_loops = dictLoops.ix[goodLoops]

        (M, C, m_count, c_count) = css.averageLoopsWithControl(these_loops, cooler_file, numShifted = 100, pad = 10)

        stored_maps['{0}-{1}bp'.format(lo, hi)] = (file_name, 'loops', M, C, m_count, c_count)
    
    return stored_maps


if __name__ == "__main__":
    main()

    
