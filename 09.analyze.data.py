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

OUTPUT_EXT = '.{resolution}{output_ext}'
TAG = 'tag1'

COOL_EXT = '.cool'

SADDLE_DIM = 5
CMD = '{snhic_path}/bin/analyze.snhic.sample.py --resolution {resolution} --genome-path {genome_path} --genome-file {genome_file} --loops-file {loops_file} --domains-file {domains_file} --scaling-ref-file {scaling_ref_file} {balance} {saddle_plots} {overwrite} {cool_file} {save_folder} {chrom_sizes_file}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--samples-file', 
               type = str,
               required = True)
@click.option('--output-dir',
              type = str, 
              required = True)
@click.option('--chrom-sizes-file',
              type = str, 
              required = True)
@click.option('--tags', 
               required = True,
               type = str, 
               help = 'The tag names seperated by a comma for the grouping of the samples.')
@click.option('--cooler-dir', 
               required = True,
               type = str, 
               help = 'The directory with the merged cooler files.')
@click.option('--cooler-ext', 
               default = COOL_EXT,
               show_default = True,
               type = str, 
               help = 'The cooler file extension.')
@click.option('--scaling-ref-file', 
               required = True,
               type = str, 
               help = 'A merged cool file (.mcool) with the scaling that can be used as a reference.')
@click.option('--resolution', 
               default = 10000, 
               type = int, 
               help = 'The resolution of the .cool files that will be analzed.')
@click.option('--genome-path', 
               required = True, 
               type = str, 
               help = 'The directory with the chomosome sequences in fasta format.')
@click.option('--genome-file', 
               required = True, 
               type = str, 
               help = 'The file with the chomosome sequences in fasta format.')
@click.option('--loops-file', 
               required = True, 
               type = str, 
               help = 'The tab delimitted file with the genomic loop positions.')
@click.option('--domains-file', 
               required = True, 
               type = str, 
               help = 'The tab delimitted file with the genomic topologically associating domain (TAD) positions.')
@click.option('--snhic-path',  
               type = str,
               default = None,
               help = 'The path to the snHi-C pipeline resources. If $SNHIC_PATH is set this is not necessary.')
@click.option('--balance/--no-balance',  
               is_flag = True,
               help = 'Set if the contact matrices were balanced. Default: --no-balance')
@click.option('--saddle-plots/--no-saddle-plots',  
               is_flag = True,
               help = 'Create the saddle plots. Default: --no-saddle-plots')
@click.option('--overwrite',  
               is_flag = True,
               help = 'Set if the previously generated files should be overwritten. Default: not set')
def main(samples_file, output_dir, chrom_sizes_file, tags, cooler_dir, cooler_ext, scaling_ref_file, resolution, genome_path, genome_file, loops_file, domains_file, snhic_path, balance, saddle_plots, overwrite):
    '''Prepare the data for future use in generating the figures.

    EXISTING_CSV: A comma seperated values (csv) file with the cellIDs of the Hi-C experiments and their metadata.

    SAVE_FOLDER: the directory that contains all the .cool files with the given resolution that will be analyzed.

    CHROM_SIZES_FILE: Chromosome order used to flip interchromosomal mates: path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose first column lists scaffol names. Any scaffolds not listed will be ordered lexicographically following the names provided.
    
    The EXISTING_CSV file is generated by 07_prepare_tables.py

    Both --genome-path and --genome-file are necessary due to different packages being used. This should be improved in future versions.
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
    for key in grouped_dict.keys():
        if isinstance(key, str) or isinstance(key, int):
            cooler_file_name = '{}{}'.format(key, OUTPUT_EXT.format(resolution = resolution, output_ext = cooler_ext)) 
        else:
            arr = list()
            for k in key:
                if not k == '':
                    arr.append(k)
                    
            cooler_file_name = '{}{}'.format('_'.join([str(e) for e in arr]), OUTPUT_EXT.format(resolution = resolution, output_ext = cooler_ext))
            
        cooler_file = os.path.join(cooler_dir, cooler_file_name)
        if not os.path.exists(cooler_file):
            sys.stderr.write('warning: merged cooler missing {}\n'.format(cooler_file))
            continue
        
        overwrite_flag = ''
        if overwrite:
            overwrite_flag = '--overwrite'

        balance_flag = ''
        if balance:
            balance_flag = '--balance'

        saddle_plots_flag = ''
        if saddle_plots:
            saddle_plots_flag = '--saddle-plots'

        cmd = CMD.format(resolution = resolution, snhic_path = snhic_path, genome_path = genome_path, genome_file = genome_file, loops_file = loops_file, domains_file = domains_file, scaling_ref_file = scaling_ref_file, balance = balance_flag, saddle_plots = saddle_plots_flag, overwrite = overwrite_flag, cool_file = cooler_file, save_folder = output_dir, chrom_sizes_file = chrom_sizes_file)
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