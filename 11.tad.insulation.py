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

TAG = 'tag1'
DISTANCE = 200000
DIAGONAL = 1
DIAMOND_SIZE = 40000
COOLER_EXT = '.cool'
OUTPUT_EXT = '.tad.insul.txt'
COOLER_FILE = '{group}.*{cooler_ext}'


CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--cooler-dir', 
                required = True, 
                type = click.Path(exists = True, resolve_path = True),
                help = 'The directory with the cooler files.')
@click.option('--domains-file', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The file with the TAD coordinates.')
@click.option('--samples-file', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The CSV file with the samples.')
@click.option('--output-dir', 
                required = True, 
                type = click.Path(dir_okay = True, resolve_path = True),
                help = 'The output directory.')
@click.option('--tags', 
               default = TAG,
               type = str, 
               help = 'Comma seperated lists of the tags used to group the samples. Default: {}'.format(TAG))
@click.option('--distance', 
                default = DISTANCE,
                type = int,
                help = 'The distance of the area around each TAD position. Default: {}'.format(DISTANCE))
@click.option('--diamond-size', 
                default = DIAMOND_SIZE,
                type = int,
                help = 'The size of the matrix for the insulation positions. Default: {}'.format(DIAMOND_SIZE))
@click.option('--diagonals', 
                default = DIAGONAL,
                type = click.IntRange(min = 0),
                help = 'The amount of diagonals that are masked in the reference cooler, where 1 is the main diagonal. Default: {}'.format(DIAGONAL))
@click.option('--cooler-ext', 
               default = COOLER_EXT,
               type = str, 
               help = 'The cooler file extension. Default: {}'.format(COOLER_EXT))
               
def main(cooler_dir, domains_file, samples_file, output_dir, tags, distance, diamond_size, diagonals, cooler_ext):
    #print(csv_file)
    df = read_samples_file(samples_file)
    if len(df) == 0:
        sys.stderr.write('error: --samples-file argument of the incorrect file type. Only support for .csv, .xls, .xlsx.')
        sys.exit()

    df = df.fillna('')

    # pool data based on sorted tags (TWO LEVEL)
    tag_names = tags.split(',')
    df = df.sort_values(tag_names)
    grouped_dict = df.groupby(tag_names).groups
    #print(grouped_dict)
    
    cooler_files = get_cooler_files(grouped_dict, cooler_dir, cooler_ext)
    #print(cooler_files)
    domain_positions = parse_domains(domains_file)

    for cooler_file in cooler_files:
        (_, cooler_file_name) = os.path.split(cooler_file)
        #print(cooler_file_name)
        #(tad_insul_contacts, domain_pos_idx, mid_idx, res) = get_tad_insulation(cooler_file, domain_positions, distance, diamond_size, diagonals)
        (tad_insul_contacts, res) = get_tad_insulation(cooler_file, domain_positions, distance, diamond_size, diagonals)
        
        #tad_insul_contacts = np.array(tad_insul_contacts)
        print(tad_insul_contacts.shape)
        #s = sum(tad_insul_contacts.transpose())
        s = sum(tad_insul_contacts)
        print(s)
        a = np.mean(s)
        stdev = np.std(tad_insul_contacts.transpose())


        idx_mid = int(len(s) / 2)

        insul = (s / min(s[idx_mid - 2 : idx_mid + 3])) - 1

        #domain_pos_idx = np.index(s[idx_mid - 2 : idx_mid + 3])
        domain_pos_idx = [i for (i, v) in enumerate(s) if v == min(s[idx_mid - 2 : idx_mid + 3])][0]
        print("domain_pos_idx: {}".format(domain_pos_idx))

        arr = list()
        arr.append('position_{}\tinsul_val_{}'.format(cooler_file_name, cooler_file_name))
        for (idx, insul_value) in enumerate(insul):
            position = (idx - domain_pos_idx) * res
            arr.append('{}\t{}'.format(position, insul_value))

        output_file_name = cooler_file_name.replace(COOLER_EXT, OUTPUT_EXT)
        output_file = os.path.join(output_dir, output_file_name)
        
        f = open(output_file, 'w')
        f.write('\n'.join(arr))
        f.close()

        print(output_file)

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


def get_tad_insulation(cooler_file, domain_positions, distance, diamond_size, diag_mask_size):
    c = cooler.Cooler(cooler_file)
    bin_size = c.info['bin-size']
    if bin_size > diamond_size:
        sys.stdout.write('warning: diamond-size ({}) changed to larger bin-size ({})\n'.format(diamond_size, bin_size))
        diamond_size = bin_size

    # if the diamond-size is not a mutiple of the bin-size change it
    if not diamond_size % bin_size == 0:
        bins_per_diamond = int(diamond_size / bin_size) + 1
        sys.stdout.write('warning: diamond-size ({}) changed to multiple of bin-size ({})\n'.format(diamond_size, bins_per_diamond))
    else:
        bins_per_diamond = int(diamond_size / bin_size)

    # if the distance is not a mutiple of the bin-size change it
    if not distance % bin_size == 0:
        bins_per_distance = int(distance / bin_size) + 1
        sys.stdout.write('warning: distance ({}) changed to multiple of bin-size ({})\n'.format(distance, bins_per_distance))
    else:
        bins_per_distance = int(distance / bin_size)

    domain_count = _count_domain_positions(domain_positions, c.chromnames)
    contacts = np.zeros((domain_count, 2 * bins_per_distance + 1))

    idx_domain = -1
    for (chrom, positions) in domain_positions.items():
        
        if chrom == 'chrM':
            continue

        if not chrom in c.chromnames:
            continue

        idx_domain += 1

        M = c.matrix(balance = False, sparse = True).fetch(chrom).toarray()
        # mask the diagonals in the contact matrix
        for diag in range(diag_mask_size):
            if diag == 0:
                np.fill_diagonal(M, 0)
            else:
                np.fill_diagonal(M[diag:], 0)
                np.fill_diagonal(M[:,diag:], 0)
        
        for pos in positions:
            mid_idx = int(pos / bin_size)  - 1 # comment out the -1 to shift one to the right
            min_idx = mid_idx - bins_per_distance
            max_idx = mid_idx + bins_per_distance

            """
            print(min_idx)
            print(mid_idx)
            print(max_idx)
            print()
            """

            idx_score = 0
            for idx in range(int(min_idx), int(max_idx) + 1):
                if idx < 0:
                    idx_score += 1
                    continue

                if idx >= len(M):
                    idx_score += 1
                    continue

                _M = get_sub_matrix(M, idx, bins_per_diamond)
                #print(_M)
                contact_count = sum(sum(_M))
                if contact_count <= 0:
                    idx_score += 1
                    continue

                contacts[idx_domain][idx_score] = contact_count
                idx_score += 1

    return (contacts, bin_size)
    #return (contacts, domain_pos_idx, mid_idx, bin_size)

def _count_domain_positions(domain_positions, chromnames):
    count = 0
    for (chrom, positions) in domain_positions.items():
        if chrom == 'chrM':
                continue

        if not chrom in chromnames:
            continue
    
        count += len(positions)

    return count


def get_sub_matrix(M, diag_pos, bins_per_diamond):
    # the n x n matrix above (i, i) of M 
    # [i + 1 - n : i + 1, i : i + n]
    # i := diag_pos
    # n := bins_per_diamond
    return M[max(0, diag_pos + 1 - bins_per_diamond) : diag_pos + 1, diag_pos : min(len(M)- 1, diag_pos + bins_per_diamond)]

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

def parse_domains(domains_file):
    domains_df = pd.read_csv(domains_file, sep = '\t')
    
    domains_df = domains_df.sort_values('chr1') # sort by the "chr1" column of the domains file
    grouped_dict = domains_df.groupby('chr1').groups

    positions = dict()
    for chrom in grouped_dict.keys():
        x1 = domains_df.ix[grouped_dict[chrom]].x1.values
        x2 = domains_df.ix[grouped_dict[chrom]].x2.values
        pos = list()
        pos.extend(list(x1))
        pos.extend(list(x2))
        pos.sort()

        positions[chrom] = pos

    #for (chrom, data) in positions.items():
    #    print('{}: {}'.format(chrom, len(data)))

    return positions
        
        

if __name__ == '__main__':
    main()