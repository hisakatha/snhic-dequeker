#! /usr/bin/env python

__version__ = '0.0.1-dev'

import os
import sys
import glob
import click

CMD = "cat {input_file} | sed '/^#/!s/#/N/g' | cooler cload pairs --metadata {stats_file} -c1 2 -p1 3 -c2 4 -p2 5 {chrom_sizes_file}:{resolution} - {cooler_file}"

INPUT_EXT = '.lowcov.pairs'
STATS_EXT = '.stats.json'
COOLER_EXT = '.cool'
COOLER_CLOAD_EXT = '.{resolution}.{cooler_ext}'

RESOLUTIONS = '1000,5000,10000,50000,100000,500000,1000000'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.argument('input_folder',
              type = str, 
              required = True)
@click.argument('output_folder',
              type = str, 
              required = True)
@click.argument('stats_folder',
              type = str, 
              required = True)
@click.argument('chrom_sizes_file',
              type = str, 
              required = True)
@click.option('--resolutions', 
               type = str, 
               default = RESOLUTIONS,
               help = 'The resolutions of the cooler files. Default: {}'.format(RESOLUTIONS))
@click.option('--input-ext', 
               default = INPUT_EXT, 
               type = str, 
               help = 'The file extension of the input pairsam files. Default: {}'.format(INPUT_EXT))
@click.option('--stats-ext', 
               default = STATS_EXT, 
               type = str, 
               help = 'The file extension of the json file with the pairsamtools statistics. Default: {}'.format(STATS_EXT))
@click.option('--output-ext', 
               default = COOLER_EXT, 
               type = str, 
               help = 'The file extension of the output cooler files. Default: {}'.format(COOLER_EXT))
def main(input_folder, output_folder, stats_folder, chrom_sizes_file, resolutions, input_ext, stats_ext, output_ext):
    '''Create the cooler cload command for the snHiC pipeline.

    INPUT_FOLDER: The directory with the pairsam files.

    OUTPUT_FOLDER: The folder with the resulting cooler files.

    STATS_FOLDER: The folder with the json file with the pairsamtools statistics.
    
    CHROM_SIZES_FILE: Chromosome order used to flip interchromosomal mates: path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose first column lists scaffol names. Any scaffolds not listed will be ordered lexicographically following the names provided.

    '''
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    input_ext = INPUT_EXT.format(input_ext)
    input_files = glob.glob(os.path.join(input_folder, '*{}'.format(input_ext)))


    cmds = list()
    for input_file in input_files:
        for resolution in resolutions.split(','):
            (_, input_file_name) = os.path.split(input_file)

            cooler_file = os.path.join(output_folder, input_file_name.replace(input_ext, '.{resolution}{output_ext}'.format(resolution = resolution, output_ext = output_ext)))
            if os.path.exists(cooler_file):
                continue

            stats_file = os.path.join(stats_folder, input_file_name.replace(input_ext, stats_ext))

            cmd = CMD.format(input_file = input_file, stats_file = stats_file, chrom_sizes_file = chrom_sizes_file, resolution = resolution, cooler_file = cooler_file)
            cmds.append(cmd)

    print('\n'.join(cmds))
    sys.stderr.write('{} commands created\n'.format(len(cmds)))


if __name__ == "__main__":
    main()
