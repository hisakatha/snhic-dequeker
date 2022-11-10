#! /usr/bin/env python
# (c) 2019. All Rights Reserved.
# 
# Code written by: Anton Goloborodko from https://github.com/mirnylab/cooltools/blob/master/doc/examples/pileups-example.ipynb
# Code altered by: Sean Powell (sean.powell@imba.oeaw.ac.at)
#

__version__ = '0.0.1-dev'

import os
import sys
import glob
import click

import cooler
import bioframe
import pybedtools
import multiprocess
import numpy as np
import pandas as pd

import cooltools
import cooltools.expected
import cooltools.snipping

import matplotlib
import matplotlib.gridspec
import matplotlib.pyplot as plt

TAG = 'tag1'
THREADS = 1
RESOLUTION = 10000
HUMAN_CHROMS = 23
MOUSE_CHROMS = 20

MCOOL_FILE = '{group}.*{mcool_ext}'
MCOOL_EXT = '.mcool'

OUTPUT_FILE_NAME = '{figure_name}.pdf'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--mcool-dir', 
                required = True, 
                type = click.Path(exists = True, resolve_path = True),
                help = 'The directory with the mcool files.')
@click.option('--samples-file', 
                required = True, 
                type = click.Path(exists = True, readable = True),
                help = "The file with experiment's the samples.")
@click.option('--tags', 
               default = TAG,
               show_default = True,
               type = str, 
               help = 'Comma seperated lists of the tags used to group the samples.')
@click.option('--figure-name', 
               required = True,
               type = str, 
               help = 'The name of the figure.')
@click.option('--genome-name', 
               required = True,
               type = click.Choice(['hg19', 'hg38', 'mm9', 'mm10']),
               help = 'The name of the reference genome.')
@click.option('--resolution', 
               default = RESOLUTION,
               show_default = True,
               type = str, 
               help = 'The resolution to be used.')
@click.option('--mcool-ext', 
               default = MCOOL_EXT,
               show_default = True,
               type = str, 
               help = 'The mcool file extension.')
@click.option('--peaks-bed', 
               required = True,
               type = click.Path(exists = True, readable = True), 
               help = 'The BED file with the ChIP-seq peaks.')
@click.option('--motifs-bed', 
               required = True,
               type = click.Path(exists = True, readable = True), 
               help = 'The BED file with the ChIP-seq motifs.')
@click.option('--output-dir', 
               required = True,
               type = click.Path(exists = True, readable = True), 
               help = 'The output path.')
@click.option('--threads', 
               default = THREADS,
               show_default = True,
               type = int, 
               help = 'The amount of threads to be used.')
def main(mcool_dir, samples_file, tags, figure_name, genome_name, resolution, mcool_ext, peaks_bed, motifs_bed, output_dir, threads):
    plt.rcParams['pdf.fonttype'] = 'truetype'
    plt.rcParams['svg.fonttype'] = 'none' # No text as paths. Assume font installed.

    plt.rcParams['font.serif'] = ['Times New Roman']
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['font.size'] = 14
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['text.usetex'] = False

    chromsizes = bioframe.fetch_chromsizes(genome_name)
    cens = bioframe.fetch_centromeres(genome_name)
    cens.set_index('chrom', inplace=True)
    cens = cens.mid

    if genome_name in ['hg19', 'hg38']:
        good_chroms = list(chromsizes.index[:HUMAN_CHROMS])

    elif genome_name in ['mm9', 'mm10']:
        good_chroms = list(chromsizes.index[:MOUSE_CHROMS])

    arms = [arm for chrom in good_chroms for arm in ((chrom, 0, cens.get(chrom,0)), (chrom, cens.get(chrom,0), chromsizes.get(chrom,0)))]
    arms = pd.DataFrame(arms, columns=['chrom','start', 'end'])

    mcool_files = get_cooler_files(samples_file, tags, mcool_dir, mcool_ext)

    for mcool_file in mcool_files:
        #MCOOL = mcool_file
        #BINSIZ   E = resolution
        clr = cooler.Cooler(f'{mcool_file}::/resolutions/{resolution}')
        with multiprocess.Pool(threads) as pool:
            expected = cooltools.expected.diagsum(clr, 
                                          list(arms.itertuples(index = False, name  =None)), 
                                          #transforms={ 'balanced': lambda p: p['count'] * p['weight1'] * p['weight2']}, 
                                          map = pool.map)

            expected_df = pd.concat([exp.reset_index().assign(chrom=reg[0], start=reg[1], end=reg[2]) for (reg, exp) in expected.items()])

            expected_df = expected_df.groupby(('chrom','diag')).aggregate({
                                                                    'n_valid':'sum',
                                                                    'count.sum':'sum',
                                                                    #'balanced.sum':'sum'
                                                                    }).reset_index()

            expected_df['count.avg'] = expected_df['count.sum'] / expected_df['n_valid']
            ctcf_peaks = pybedtools.BedTool(peaks_bed).sort()
            #########
            #ctcf_motifs = pybedtools.BedTool('./encode_motifs.hg38.ctcf_known1.liftover.bed.gz').sort()
            #########
            ctcf_motifs = pybedtools.BedTool(motifs_bed).sort()
            ctcf_motifs_w_peaks = ctcf_motifs.intersect(ctcf_peaks).to_dataframe()
            ctcf_motifs_w_peaks['mid'] = (ctcf_motifs_w_peaks.start + ctcf_motifs_w_peaks.end) / 2\

            WINDOW_HALF_SIZE = 100000
            snipping_windows = cooltools.snipping.make_bin_aligned_windows(resolution, ctcf_motifs_w_peaks.chrom.values, ctcf_motifs_w_peaks.mid.values, WINDOW_HALF_SIZE)
            oe_snipper = cooltools.snipping.ObsExpSnipper(clr, expected_df)

            with multiprocess.Pool(threads) as pool:
                oe_pile = cooltools.snipping.pileup(snipping_windows, oe_snipper.select, oe_snipper.snip, map = pool.map)

            collapsed_pile_plus = np.nanmean(oe_pile[:, :, ctcf_motifs_w_peaks.strand=='+'], axis = 2)
            collapsed_pile_minus = np.nanmean(oe_pile[:, :, ctcf_motifs_w_peaks.strand=='-'], axis = 2)
            
            plt.imshow(np.log2(collapsed_pile_plus), vmax = 0.5, vmin = -0.5, cmap='coolwarm')
            plt.colorbar(label = 'log2 mean obs/exp')

            ticks_pixels = np.linspace(0, WINDOW_HALF_SIZE * 2 // resolution, 5)
            ticks_kbp = ((ticks - ticks[-1] / 2) * resolution // 1000).astype(int)
            
            plt.xticks(ticks_pixels, ticks_kbp)
            plt.yticks(ticks_pixels, ticks_kbp)
            plt.xlabel('relative position, kbp')
            plt.ylabel('relative position, kbp')


def get_cooler_files(samples_file, tags, cooler_dir, cooler_ext):
    tag_names = tags.split(',')

    df = read_samples_file(samples_file)
    df = df[df['include'] == True]
    if len(df) == 0:
        sys.stderr.write('error: --samples-file argument of the incorrect file type. Only support for .csv, .xls, .xlsx.')
        sys.exit()

    df = df.fillna('')
    df = df.sort_values(tag_names)
    grouped_dict = df.groupby(tag_names).groups

    file_names = list()
    for key in grouped_dict.keys():
        if isinstance(key, str) or isinstance(key, int):
            cooler_file_name = MCOOL_FILE.format(group = key, mcool_ext = cooler_ext)
        else:
            arr = list()
            for k in key:
                if not k == '':
                    arr.append(k)

            cooler_file_name = MCOOL_FILE.format(group = '_'.join(arr), mcool_ext = cooler_ext)
        
        cooler_files = glob.glob(os.path.join(cooler_dir, cooler_file_name))
        if len(cooler_files) == 0:
            continue
    
        file_names.append(cooler_files[0])
        
    return file_names

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