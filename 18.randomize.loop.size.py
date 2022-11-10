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

OUTPUT_FILE_NAME = '{condition_a}.{condition_b}.loop.strength.random.txt'
COOLER_EXT = '.{resolution}{cooler_ext}'
TAG = 'tag1'
COOLER_EXT = '.cool'
SMOOTH_SIGMA = 2
MAX_EXP_LOOP_SIZE = 1e6
#####
# needs to be changed for the paired bootstrap analysis
CMD = '{snhic_path}/bin/randomize.loop.size.py --condition-a {condition_a} --condition-b {condition_b} --sample-size-a {sample_size_a} --sample-size-b {sample_size_b} --sample-ids {sample_ids} --output-file {output_file} --genome-path {genome_path} --cooler-dir {cooler_dir} --resolution {resolution} --threads {threads} --iterations {iterations}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--samples-file', 
               type = str,
               required = True)
@click.option('--output-dir',
              type = str, 
              required = True)
@click.option('--tags', 
               required = True,
               type = str, 
               help = 'The tag names seperated by a comma for the grouping of the samples. This should correspond to the conditions.')
@click.option('--conditions', 
               required = True,
               type = str, 
               help = 'The condition pairs. The pairs are seperated by a semicolon (;) while the conditions are seperated by a comma (,). e.g. cond1,cond2;cond3,cond4')
@click.option('--cooler-dir', 
               required = True,
               type = str, 
               help = 'The directory with the merged cooler files.')
@click.option('--cooler-ext', 
               default = COOLER_EXT,
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
               help = 'The maximun loop size taken into account.')
@click.option('--snhic-path',  
               type = str,
               default = None,
               help = 'The path to the snHi-C pipeline resources. If $SNHIC_PATH is set this is not necessary.')
def main(samples_file, output_dir, tags, conditions, cooler_dir, cooler_ext, resolution, threads, iterations, genome_path, ref_cooler, chroms_file, smooth_sigma, max_loop_size, snhic_path):
    '''Create randomization command for the aggregate loop analysis.

    '''

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    assert(os.path.exists(output_dir))

    tag_names = tags.split(',')
    df = read_samples_file(samples_file)
    df = df[df['include'] == True]

    if len(df) == 0:
        sys.stderr.write('error: --samples-file argument of the incorrect file type. Only support for .csv, .xls, .xlsx.')
        sys.exit()

    df = df.fillna('')
    df = df.sort_values(tag_names)

    if snhic_path == None:
        snhic_path = os.environ['SNHIC_PATH']

    grouped_dict = df.groupby(tag_names).groups
    cmds = list()
    
    condition_pairs = conditions.split(';')
    for condition_pair in condition_pairs:
        (condition_a, condition_b) = condition_pair.split(',')

        cellIDs_a = df.ix[grouped_dict[condition_a]].cellID.values
        cellIDs_b = df.ix[grouped_dict[condition_b]].cellID.values
        cellIDs = list()
        cellIDs.extend(cellIDs_a)
        cellIDs.extend(cellIDs_b)

        output_file_name = OUTPUT_FILE_NAME.format(condition_a = condition_a, condition_b = condition_b)
        output_file = os.path.join(output_dir, output_file_name)

        cmd = CMD.format(snhic_path = snhic_path, condition_a = condition_a, condition_b = condition_b, sample_size_a = len(cellIDs_a), sample_size_b = len(cellIDs_b), sample_ids = ','.join([str(i) for i in cellIDs]), output_file = output_file, genome_path = genome_path, cooler_dir = cooler_dir, resolution = resolution, threads = threads, iterations = iterations)
        cmds.append(cmd)

    print('\n'.join(cmds))
    sys.stderr.write('{} commands created\n'.format(len(cmds)))

def read_samples_file(samples_file):
    (_, file_ext) = os.path.splitext(samples_file)

    if file_ext == '.csv':
        df = pd.read_csv(samples_file)
    elif file_ext == '.xls':
        df = pd.read_excel(samples_file)
    elif file_ext == '.xlsx':
        df = pd.read_excel(samples_file)
    else:
        df = None

    return df

if __name__ == '__main__':
    main()