#! /usr/bin/env python
# (c) 2018. All Rights Reserved.
# Code written by: Sean Powell (sean.powell@imba.oeaw.ac.at)

__version__ = '0.0.1-dev'

import os
import click
import cooler
import pickle
import pandas as pd
from snHiC.lib import coolerAnalysisShared as css
from snHiC.lib import analyzeSnHiCSample as ass

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'], }
@click.version_option(version=__version__)
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('cool_file',
                type=str,
                required=True)
@click.argument('save_folder',
                type=str,
                required=True)
@click.argument('chrom_sizes_file',
                type=str,
                required=True)
@click.option('--scaling-ref-file',
              type=str,
              help='A cool file (.cool) to be used as a scaling reference.')
@click.option('--resolution',
              default=10000,
              type=int,
              help='The resolution of the .cool files that will be analzed.')
@click.option('--genome-path',
              required=True,
              type=str,
              help='The directory with the chomosome sequences in fasta format.')
@click.option('--genome-file',
              required=True,
              type=str,
              help='The file with the chomosome sequences in fasta format.')
@click.option('--loops-file',
              required=True,
              type=str,
              help='The tab delimitted file with the genomic loop positions.')
@click.option('--domains-file',
              required=True,
              type=str,
              help='The tab delimitted file with the genomic topologically associating domain (TAD) positions.')
@click.option('--balance',
              is_flag=True,
              help='Set if the contact matrices were balanced. Default: not set')
@click.option('--saddle-plots',
              is_flag=True,
              help='Create the saddle plots. Default: not set')
@click.option('--overwrite',
              is_flag=True,
              help='Set if the previously generated files should be overwritten. Default: not set')
# make plots for figure 1
def main(cool_file, save_folder, chrom_sizes_file, scaling_ref_file, resolution, genome_path, genome_file, loops_file, domains_file, balance, saddle_plots, overwrite):
    '''Prepare the data for future use in generating the figures.

    COOL_FILE: The .cool file of the sample to be analyzed.

    SAVE_FOLDER: the directory that contains all the .cool files with the given resolution that will be analyzed.

    CHROM_SIZES_FILE: Chromosome order used to flip interchromosomal mates: path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose first column lists scaffol names. Any scaffolds not listed will be ordered lexicographically following the names provided.

    '''

    assert(os.path.exists(cool_file))
    assert(os.path.exists(genome_path))
    assert(os.path.exists(loops_file))
    assert(os.path.exists(domains_file))
    assert(os.path.exists(chrom_sizes_file))

    loopAndDomainDict = css.getLoopsAndDomains(genome_path, loops_file, domains_file)
    ass.deploySample(cool_file, save_folder, genome_path, scaling_ref_file, loopAndDomainDict, saddle_plots, overwrite)


if __name__ == '__main__':
    main()
