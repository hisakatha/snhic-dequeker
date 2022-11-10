#! /usr/bin/env python
# (c) 2019. All Rights Reserved.
# Code written by: 
#       Maxim
#       Hugo
#       Sean Powell (sean.powell@imba.oeaw.ac.at)

__version__ = '0.0.1-dev'

import os
import sys
import click
import pandas as pd

OUTPUT_FILE_NAME = '{prefix}.loop.strength.bs.txt'
COOLER_FILE = '{condition}.{resolution}{cooler_ext}'
TAG = 'tag1'
COOL_EXT = '.cool'
#####
# needs to be changed for the paired bootstrap analysis
CMD = '{snhic_path}/bin/bootstrap.loop.strength.py --condition-a {condition_a} --condition-b {condition_b} --cooler-file-a {cooler_file_a} --cooler-file-b {cooler_file_b} --output-file {output_file} --genome-path {genome_path} --loops-file {loops_file} --threads {threads} --iterations {iterations}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--output-dir',
              type = str, 
              required = True)
@click.option('--conditions', 
               required = True,
               type = str, 
               help = 'The condition pairs. The pairs are seperated by a semicolon (;) while the conditions are seperated by a comma (,). e.g. cond1,cond2;cond3,cond4')
@click.option('--cooler-dir', 
               required = True,
               type = str, 
               help = 'The directory with the merged cooler files.')
@click.option('--cooler-ext', 
               default = COOL_EXT,
               show_default = True,
               type = str, 
               help = 'The cooler file extension.')
@click.option('--resolution', 
               default = 10000,
               show_default = True,
               type = int, 
               help = 'The resolution of the .cool files that will be analzed.')
@click.option('--threads', 
               default = 8, 
               show_default = True,
               type = int, 
               help = 'The concurrent threads that will be started for the bootstrapping.')
@click.option('--iterations', 
               default = 100, 
               show_default = True,
               type = int, 
               help = 'The iterations that will be bootstrapped.')
@click.option('--genome-path', 
               required = True, 
               type = str, 
               help = 'The directory with the chomosome sequences in fasta format.')
@click.option('--loops-file', 
               required = True, 
               type = str, 
               help = 'The tab delimitted file with the genomic loop positions.')
@click.option('--snhic-path',  
               type = str,
               default = None,
               help = 'The path to the snHi-C pipeline resources. If $SNHIC_PATH is set this is not necessary.')
def main(output_dir, conditions, cooler_dir, cooler_ext, resolution, threads, iterations, genome_path, loops_file, snhic_path):
    '''Create bootstrapping command for the aggregate loop analysis.

    '''
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    assert(os.path.exists(output_dir))

    if snhic_path == None:
        snhic_path = os.environ['SNHIC_PATH']

    cmds = list()
    condition_pairs = conditions.split(';')
    for condition_pair in condition_pairs:
        (condition_a, condition_b) = condition_pair.split(',')
        cooler_file_name_a = COOLER_FILE.format(condition = condition_a, resolution = resolution, cooler_ext = cooler_ext)
        cooler_file_a = os.path.join(cooler_dir, cooler_file_name_a)
        cooler_file_name_b = COOLER_FILE.format(condition = condition_b, resolution = resolution, cooler_ext = cooler_ext)
        cooler_file_b = os.path.join(cooler_dir, cooler_file_name_b)

        output_prefix = '{}.{}'.format(condition_a, condition_b)
        output_file_name = OUTPUT_FILE_NAME.format(prefix = output_prefix)
        output_file = os.path.join(output_dir, output_file_name)

        cmd = CMD.format(snhic_path = snhic_path, condition_a = condition_a, condition_b = condition_b, cooler_file_a = cooler_file_a, cooler_file_b = cooler_file_b, output_file = output_file, genome_path = genome_path, loops_file = loops_file, threads = threads, iterations = iterations)
        cmds.append(cmd)

    print('\n'.join(cmds))
    sys.stderr.write('{} commands created\n'.format(len(cmds)))

if __name__ == '__main__':
    main()