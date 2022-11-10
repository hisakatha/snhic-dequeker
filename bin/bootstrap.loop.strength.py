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
import click

import pandas as pd
import numpy as np
import multiprocessing as mp

from scipy import stats

from pandas import ExcelWriter
from mirnylib.genome import Genome
from snHiC.lib import coolerAnalysisShared as css


CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--condition-a', 
               type = str,
               required = True)
@click.option('--condition-b', 
               type = str,
               required = True)
@click.option('--cooler-file-a', 
               type = click.Path(exists = True, resolve_path = True),
               required = True)
@click.option('--cooler-file-b', 
               type = click.Path(exists = True, resolve_path = True),
               required = True)
@click.option('--output-file',
              type = str, 
              required = True,
              help = 'The name of the outfile in .xlsx format.')
@click.option('--genome-path', 
               required = True, 
               type = click.Path(exists = True, resolve_path = True), 
               help = 'The directory with the chomosome sequences in fasta format.')
@click.option('--loops-file', 
               required = True, 
               type = click.Path(exists = True, resolve_path = True), 
               help = 'The tab delimitted file with the genomic loop positions.')
@click.option('--threads', 
               default = 8, 
               show_default = True,
               type = int, 
               help = 'The concurrent threads that will be started for the bootstrapping.')
@click.option('--iterations', 
               default = 100, 
               show_default = True,
               type = int, 
               help = 'The iterations of the bootstrapping.')
def main(condition_a, condition_b, cooler_file_a, cooler_file_b, output_file, genome_path, loops_file, threads, iterations):
    '''Bootstrap the aggregate loops and domains analysis.
    '''
    gen = Genome(genome_path, readChrms = ["#", "X"])
    loops_df = pd.read_csv(loops_file, sep = "\t")

    if threads <= 0:
        threads = 1

    if threads > iterations:
        threads = iterations

    print(os.path.split(cooler_file_a)[1])
    print(os.path.split(cooler_file_b)[1])

    results = get_results(loops_df, gen, cooler_file_a, cooler_file_b, threads, iterations)
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

def get_results(loops_df, gen, cooler_file_a, cooler_file_b, max_concurrent, jobs):
    manager = mp.Manager()
    return_list = manager.list()
    running_jobs = list()
    while True:
        if len(return_list) < jobs:
            for i in range(max_concurrent - len(running_jobs)):
                p = mp.Process(target = getAggrLoopStrength, args = (loops_df, gen, cooler_file_a, cooler_file_b, return_list))
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

def getAggrLoopStrength(loops_df, gen, cooler_file_a, cooler_file_b, return_list):
    # random sample with replacement, this is basicaaly the bootstrap method
    loops_df = loops_df.sample(len(loops_df), replace = True, random_state = np.random.RandomState())
    #print(loops_df)
    
    print('getLoops')
    dictLoops = css.getLoops(gen, loops_df)
    print('doPileUp')
    stored_maps_a = doPileUp(cooler_file_a, dictLoops)
    stored_maps_b = doPileUp(cooler_file_b, dictLoops)
    print('showLoopPileUp')
    (M_a, C_a, m_count_a, c_count_a) = css.showLoopPileUp(stored_maps_a)
    (M_b, C_b, m_count_b, c_count_b) = css.showLoopPileUp(stored_maps_b)
    print('getLoopStrength')
    strength_a = getLoopStrength(M_a, m_count_a, C_a, c_count_a)
    strength_b = getLoopStrength(M_b, m_count_b, C_b, c_count_b)
    

    return_list.append((strength_a, strength_b))

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

    
