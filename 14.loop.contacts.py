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

from scipy import stats

TAG = 'tag1'
COOLER_EXT = '.cool'
OUTPUT_EXT = '.loop.contacts.txt'
COOLER_FILE = '{group}.*{cooler_ext}'


CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--cooler-dir', 
                required = True, 
                type = click.Path(exists = True, resolve_path = True),
                help = 'The directory with the cooler files.')
@click.option('--loops-file', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The file with the loop coordinates.')
@click.option('--samples-file', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The .xls|.xlsx|.csv file with the samples.')
@click.option('--output-dir', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The output directory.')
@click.option('--tags', 
               default = TAG,
               type = str, 
               show_default = True,
               help = 'Comma seperated lists of the tags used to group the samples.')
@click.option('--diagonals', 
               default = 2,
               type = int, 
               show_default = True,
               help = 'The amount of diagonals that will be masked.')
@click.option('--cooler-ext', 
               default = COOLER_EXT,
               type = str, 
               show_default = True,
               help = 'The cooler file extension.')
               
def main(cooler_dir, loops_file, samples_file, output_dir, tags, diagonals, cooler_ext):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    df = read_samples_file(samples_file)
    if len(df) == 0:
        sys.stderr.write('error: --samples-file argument of the incorrect file type. Only support for .csv, .xls, .xlsx.')
        sys.exit()

    df = df.fillna('')

    # pool data based on sorted tags (TWO LEVEL)
    tag_names = tags.split(',')
    df = df.sort_values(tag_names)
    grouped_dict = df.groupby(tag_names).groups
    
    cooler_files = get_cooler_files(grouped_dict, cooler_dir, cooler_ext)
    
    loop_positions = parse_loops(loops_file)
    dic_contacts_norm = dict()
    for cooler_file in cooler_files:
        (_, cooler_file_name) = os.path.split(cooler_file)

        #print(cooler_file_name)
        chrom_contact_counts = get_loop_contacts(cooler_file, loop_positions, diagonals)
        #print(chrom_contact_counts)

        arr = list()
        arr.append('chrom\tloop_start\tloop_end\tloop_start_idx\tloop_bin_length\tcontacts_{}\tcontacts_norm_{}'.format(cooler_file_name, cooler_file_name))
        
        contacts_norm_all = list()
        for (chrom, positions) in chrom_contact_counts.items():
            for pos in positions:
                (loop_start, loop_end, loop_start_idx, loop_bin_length, contacts, contacts_norm) = pos
                arr.append('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrom, loop_start, loop_end, loop_start_idx, loop_bin_length, contacts, contacts_norm))

                contacts_norm_all.append(contacts_norm)

        dic_contacts_norm[cooler_file_name] = contacts_norm_all

        output_file_name = cooler_file_name.replace(COOLER_EXT, OUTPUT_EXT)
        output_file = os.path.join(output_dir, output_file_name)

        f = open(output_file, 'w')
        f.write('\n'.join(arr))
        f.close()
    
    #print some statistics
    equal_var = True

    print('cond_a\tcond_b\twilcoxon_p_value\tt_test_p_value')
    for cooler_file_name_a in dic_contacts_norm:
        for cooler_file_name_b in dic_contacts_norm:
            if cooler_file_name_a == cooler_file_name_b:
                continue

            #(stat, p_value) = stats.ttest_ind(dic_contacts_norm[cooler_file_name_a], dic_contacts_norm[cooler_file_name_b], equal_var = equal_var)
            (_, t_test_p_value) = stats.ttest_rel(dic_contacts_norm[cooler_file_name_a], dic_contacts_norm[cooler_file_name_b])
            (_, wilcoxon_p_value) = stats.mannwhitneyu(dic_contacts_norm[cooler_file_name_a], dic_contacts_norm[cooler_file_name_b])
            print('{}\t{}\t{}\t{}'.format(cooler_file_name_a, cooler_file_name_b, wilcoxon_p_value, t_test_p_value))


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


def get_loop_contacts(cooler_file, loop_positions, diag_mask_size):
    c = cooler.Cooler(cooler_file)
    bin_size = c.info['bin-size']
    sample_count = c.info['metadata']['files_merged']
    chrom_contact_counts = dict()

    for (chrom, positions) in loop_positions.items():
        if chrom == 'chrM':
            continue

        if not chrom in c.chromnames:
            continue
        chrom_contact_counts[chrom] = list()
        M = c.matrix(balance = False, sparse = True).fetch(chrom).toarray()

        # mask the diagonals in the contact matrix
        for diag in range(diag_mask_size):
            if diag == 0:
                np.fill_diagonal(M, 0)
            else:
                np.fill_diagonal(M[diag:], 0)
                np.fill_diagonal(M[:,diag:], 0)
        
        for pos in positions:
            min_pos = min(pos)
            max_pos = max(pos)

            min_idx = int(min_pos / bin_size) - 1
            if min_pos < bin_size:
                min_idx = 0

            max_idx = int(max_pos / bin_size) - 1
            if not max_pos % bin_size == 0:
                max_idx += 1

            idx_dist = max_idx - min_idx + 1

            _M = get_sub_matrix(M, min_idx, idx_dist)
            contact_count = sum(sum(_M))

            #print(min_pos, max_pos, min_idx, idx_dist, contact_count, contact_count / sample_count)

            chrom_contact_counts[chrom].append((min_pos, max_pos, min_idx, idx_dist, contact_count, contact_count / sample_count))


    return chrom_contact_counts


def get_sub_matrix(M, pos, length):
    # the length x length matrix from point pos 
    # [pos : pos + length, pos : pos + length]
    # 
    return M[pos : pos + length, pos : pos + length]

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
        
        cooler_files = glob.glob(os.path.join(cooler_dir, cooler_file_name))
        if len(cooler_files) == 0:
            continue
    
        file_names.append(cooler_files[0])
        
    return file_names

def parse_loops(loops_file):
    loops_df = pd.read_csv(loops_file, sep = '\t')
    
    loops_df = loops_df.sort_values('chr1') # sort by the "chr1" column of the domains file
    grouped_dict = loops_df.groupby('chr1').groups

    positions = dict()
    for chrom in grouped_dict.keys():
        x1 = loops_df.ix[grouped_dict[chrom]].x1.values
        x2 = loops_df.ix[grouped_dict[chrom]].x2.values
        y1 = loops_df.ix[grouped_dict[chrom]].y1.values
        y2 = loops_df.ix[grouped_dict[chrom]].y2.values
        pos = list(zip(x1, x2, y1, y2))
        pos.sort()

        positions[chrom] = pos

    return positions
        
        

if __name__ == '__main__':
    main()