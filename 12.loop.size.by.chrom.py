#! /usr/bin/env python
# (c) 2019. All Rights Reserved.
# Code written by: Sean Powell (sean.powell@imba.oeaw.ac.at)

__version__ = '0.0.1-dev'

import os
import sys
import glob
import click
import cooler
import numpy as np
import pandas as pd

from pandas import ExcelWriter
from snHiC.lib import coolerAnalysisShared as css

from scipy.ndimage.filters import gaussian_filter1d

TAG = 'tag1'
COOLER_EXT = '.cool'
SMOOTH_SIGMA = 2
MAX_EXP_LOOP_SIZE = 1e6
COOLER_FILE = '{group}.*{cooler_ext}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--cooler-dir', 
                required = True, 
                type = click.Path(exists = True, resolve_path = True),
                help = 'The directory with the cooler files.')
@click.option('--xlsx-file', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The .xlsx file with the samples.')
@click.option('--output-file', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The output file name as an .xlsx file.')
@click.option('--ref-cooler', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The reference cooler used to filter out noise.')
@click.option('--chroms-file', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'A file with the chromosome sizes')
@click.option('--tags', 
               default = TAG,
               type = str, 
               show_default = True,
               help = 'Comma seperated lists of the tags used to group the samples.')
@click.option('--cooler-ext', 
               default = COOLER_EXT,
               type = str, 
               show_default = True,
               help = 'The cooler file extension.')
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
@click.option('--bin-size', 
               default = None,
               type = str, 
               help = 'Instead of binning by chromosome make bin of this size. Default: bin by chromosome')        
def main(cooler_dir, xlsx_file, output_file, ref_cooler, chroms_file, tags, cooler_ext, smooth_sigma, max_loop_size, bin_size):
    #print(xlsx_file)
    #df = pd.read_csv(csv_file)
    df = pd.read_excel(xlsx_file)
    df = df.fillna('')

    # pool data based on sorted tags (TWO LEVEL)
    tag_names = tags.split(',')
    df = df.sort_values(tag_names)
    grouped_dict = df.groupby(tag_names).groups
    
    ref_c = cooler.Cooler(ref_cooler)
    ref_res = ref_c.info['bin-size']
    #print(ref_c.chromnames)

    dicChromsSize = parse_chroms_file(chroms_file)
    if dicChromsSize == None:
        sys.stderr.write('error: parsing {}\n'.format(os.path.split(chroms_file)[1]))
        sys.exit()

    cooler_files = get_cooler_files(grouped_dict, cooler_dir, cooler_ext)
    #print(cooler_files)
    column_names = list()
    data = list()
    chr_names = list()
    for cooler_file in cooler_files:
        (_, cooler_file_name) = os.path.split(cooler_file)
        print(cooler_file_name)

        column_names.append('{}_start'.format(cooler_file_name))
        column_names.append('{}_end'.format(cooler_file_name))
        column_names.append('{}_contacts'.format(cooler_file_name))
        start_loop_size = list()
        end_loop_size = list()
        chrom_contacts = list()

        c = cooler.Cooler(cooler_file)
        res = c.info['bin-size']
        if not res == ref_res:
            print('cooler resolution {} does not match the reference resolution {}, skipping...'.format(res, ref_res))
            continue

        chrom_names = css.map_chromnames(c.chromnames, ref_c.chromnames)
        if len(chr_names) == 0:
            for (chr_name, _) in chrom_names:
                chr_names.append(chr_name)

        #print(chr_names)
        for (chrom, ref_chrom) in chrom_names:
            vals_vals = get_cooler_scaling(c, ref_c, chrom, ref_chrom, dicChromsSize)
            for (pos, vals) in vals_vals.items():
                start = pos[0]
                end = pos[1]
                contacts = vals[2]
                (loop_size_start, loop_size_end) = do_slopes(vals, res, smooth_sigma, max_loop_size)

                start_loop_size.append( loop_size_start * res / 1e3)
                end_loop_size.append(loop_size_end * res / 1e3)
                chrom_contacts.append(contacts)
                #print('{}:{:.0f}-{:.0f}kbp:{}'.format(chrom, loop_size_start * res / 1e3 , loop_size_end * res / 1e3, contacts ))

        data.append(start_loop_size)
        data.append(end_loop_size)
        data.append(chrom_contacts)

    data = np.array(data)
    I = pd.Index(chr_names, name = 'chromosome')
    C = pd.Index(column_names)
    df = pd.DataFrame(data.transpose(), index = I, columns = C)

    writer = ExcelWriter(output_file)
    df.to_excel(writer)
    writer.save()

def get_cooler_scaling(c, ref_c, chrom, ref_chrom, dic_chroms_size, bin_size = None, logFactor = 1.15):
    # this function assumes the chrom is in both c and ref_c
    # therefore this needs to be checked before the function is called
    
    #Pc_list = list()
    #mids_list = list()
    if not chrom in dic_chroms_size:
        sys.stderr.write('error: cannot find {} in chromosome length file\n'.format(chrom))
        sys.exit()

    M = c.matrix(balance = False).fetch(chrom)
    M_ref = ref_c.matrix(balance = False).fetch(ref_chrom)
    
    dicScaling = dict()
    if bin_size == None:
        (_Pc, mids, contacts) = get_matrix_scaling(M, M_ref, logFactor)

        dicScaling[(1, dic_chroms_size[chrom])] = (_Pc, mids, contacts)
    else:
        pass


    return dicScaling


def get_matrix_scaling(M, M_ref, logFactor):
    marginals = np.nansum(M_ref, axis = 0)
    mask = marginals > 0
    mask2d = mask[:, None] * mask[None, :]

    (Pc, mids) = css.getMatrixScaling(M, inMask = mask2d, measureType = 'sum', logFactor = logFactor)
    #print(Pc)
    
    #Pc_list.append(Pc)
    #mids_list.append(mids)

    _Pc = np.zeros((1, len(Pc))) * np.nan
    #(si, s) = Pc
    _Pc[0, 0:len(Pc)] = Pc

    _Pc = np.nanmean(_Pc, axis = 0)
    contacts = sum(sum(M))

    return (_Pc, mids, contacts)

def do_slopes(vals, resolution, smooth_sigma, max_exp_loop_size):
    
    #print(vals)


    dlogx = np.diff(np.log(vals[1]))
    slope = np.diff(np.log(vals[0])) / dlogx
    
    loop_size_filt_idx = np.argmax(gaussian_filter1d(slope[(vals[1][1:] * resolution < max_exp_loop_size)], smooth_sigma))
    loop_size_from = vals[1][loop_size_filt_idx - 1]
    loop_size_to = vals[1][loop_size_filt_idx]

    return (loop_size_from, loop_size_to)
        

def get_cooler_files(grouped_dict, cooler_dir, cooler_ext):
    file_names = list()
    for key in grouped_dict.keys():
        if isinstance(key, str) or isinstance(key, int):
            cooler_file_name = COOLER_FILE.format(group = key, cooler_ext = cooler_ext)
        else:
            arr = list()
            for k in key:
                if not k == '':
                    arr.append(k)

            cooler_file_name = COOLER_FILE.format(group = '_'.join(arr), cooler_ext = cooler_ext)
        #print(cooler_file_name)
        
        cooler_files = glob.glob(os.path.join(cooler_dir, cooler_file_name))
        if len(cooler_files) == 0:
            continue
    
        file_names.append(cooler_files[0])
        
    return file_names

def parse_chroms_file(chroms_file):
    f = open(chroms_file)

    dicChromSize = dict()
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue

        try:
            (chrom, size) = line.split()
        except:
            continue

        try:
            size = int(size)
        except:
            return None

        dicChromSize[chrom] = size

    return dicChromSize
        
if __name__ == '__main__':
    main()
