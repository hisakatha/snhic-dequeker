#! /usr/bin/env python

__version__ = '0.0.1-dev'

import os
import sys
import glob
import click

#CMD_FILTERBYCOV = "pairtools parse --output-stats {parse_stats_file} --assembly {assembly_name} --chroms-path {chroms_path} {sam_file} | pairtools sort --nproc {threads} --memory 2G | pairtools select '(pair_type == \"UU\") or (pair_type == \"UR\") or (pair_type == \"RU\")' --output-rest >( pairtools split --output-pairs {unmapped_pairs_file} --output-sam {unmapped_sam_file} ) | pairtools dedup --output-stats {dedup_stats_file} --max-mismatch 500 --output-dups >( pairtools markasdup | pairtools split --output-pairs {dups_pairs_file} --output-sam {dups_sam_file} ) | pairtools filterbycov --output-stats {filterbycov_stats_file} --output-highcov >( pairtools split --output-pairs {highcov_pairs_file} --output-sam {highcov_sam_file} ) --output >( pairtools split --output-pairs {output_file} --output-sam {output_sam_file} )"
#CMD_DEDUP = "pairtools parse --output-stats {parse_stats_file} --assembly {assembly_name} --chroms-path {chroms_path} {sam_file} | pairtools sort --nproc {threads} --memory 2G | pairtools select '(pair_type == \"UU\") or (pair_type == \"UR\") or (pair_type == \"RU\")' --output-rest >( pairtools split --output-pairs {unmapped_pairs_file} --output-sam {unmapped_sam_file} ) | pairtools dedup --output-stats {dedup_stats_file} --max-mismatch 500 --output-dups >( pairtools markasdup | pairtools split --output-pairs {dups_pairs_file} --output-sam {dups_sam_file} ) | pairtools split --output-pairs {output_file} --output-sam {output_sam_file}"
CMD_FILTERBYCOV = "pairtools parse --output-stats {parse_stats_file} --assembly {assembly_name} --chroms-path {chroms_path} {sam_file} | pairtools sort --nproc {threads} --memory 2G | pairtools select '(pair_type == \"UU\") or (pair_type == \"UR\") or (pair_type == \"RU\")' | pairtools dedup --output-stats {dedup_stats_file} --max-mismatch 500 | pairtools filterbycov --output-stats {filterbycov_stats_file} --output >( pairtools split --output-pairs {output_file} --output-sam {output_sam_file} )"
CMD_DEDUP = "pairtools parse --output-stats {parse_stats_file} --assembly {assembly_name} --chroms-path {chroms_path} {sam_file} | pairtools sort --nproc {threads} --memory 2G | pairtools select '(pair_type == \"UU\") or (pair_type == \"UR\") or (pair_type == \"RU\")' | pairtools dedup --output-stats {dedup_stats_file} --max-mismatch 500 | pairtools split --output-pairs {output_file} --output-sam {output_sam_file}"

CHROMS_FILE = '/groups/tachibana/sean/data/genomes/mus_musculus/mm9/mm9.reduced.chrom.sizes'

PARSE_STATS_SUBFOLDER = 'stats.parse'
DEDUP_STATS_SUBFOLDER = 'stats.dedup'
FILTERBYCOV_STATS_SUBFOLDER = 'stats.filterbycov'
DEDUP_PAIRS_SUBFOLDER = 'pairs.dups'
DEDUP_SAM_SUBFOLDER = 'sam.dups'
LOWCOV_PAIRS_SUBFOLDER = 'pairs.lowcov'
LOWCOV_SAM_SUBFOLDER = 'sam.lowcov'
HIGHCOV_PAIRS_SUBFOLDER = 'pairs.highcov'
HIGHCOV_SAM_SUBFOLDER = 'sam.highcov'
UNMAPPED_PAIRS_SUBFOLDER = 'pairs.unmapped'
UNMAPPED_SAM_SUBFOLDER = 'sam.unmapped'

SAM_EXT = '.sam'
PAIRS_EXT = '.pairs'
STATS_EXT = '.stats'
#STATS_EXT = '.stats.json'

PARSE_EXT = '.parse{}'
DEDUP_EXT = '.dedup{}'
FILTERBYCOV_EXT = '.filterbycov{}'
UNMAPPED_EXT = '.unmapped{}'
DUPS_EXT = '.dups{}'
LOWCOV_EXT = '.lowcov{}'
HIGHCOV_EXT = '.highcov{}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.argument('sam_folder',
              type = str, 
              required = True)
@click.argument('output_folder',
              type = str, 
              required = True)
@click.argument('output_sam_folder',
              type = str, 
              required = True)
@click.argument('chrom_sizes_file',
              type = str, 
              required = True)
@click.option('--assembly-name', 
               type = str, 
               required = True,
               help = 'The name of the genome assembly being used. (e.g. mm9)')
@click.option('--working-dir', 
               type = str, 
               default = None,
               help = 'The working directory. Default: cwd')
@click.option('--parse-stats-folder', 
               type = str, 
               default = PARSE_STATS_SUBFOLDER,
               help = 'The name of the folder for pairsam parse stats files.')
@click.option('--dedup-stats-folder', 
               type = str, 
               default = DEDUP_STATS_SUBFOLDER,
               help = 'The name of the folder for pairsam dedup stats files.')
@click.option('--filterbycov-stats-folder', 
               type = str, 
               default = FILTERBYCOV_STATS_SUBFOLDER,
               help = 'The name of the folder for pairsam filterbycov stats files.')
@click.option('--unmapped-pairs-folder', 
               type = str, 
               default = UNMAPPED_PAIRS_SUBFOLDER,
               help = 'The name of the folder for unmapped pairsam files.')
@click.option('--unmapped-sam-folder', 
               type = str, 
               default = UNMAPPED_SAM_SUBFOLDER,
               help = 'The name of the folder for unmapped sam files.')
@click.option('--dups-pairs-folder', 
               type = str, 
               default = DEDUP_PAIRS_SUBFOLDER,
               help = 'The name of the folder for duplicates pairsam files.')
@click.option('--dups-sam-folder', 
               type = str, 
               default = DEDUP_SAM_SUBFOLDER,
               help = 'The name of the folder for duplicates sam files.')
@click.option('--lowcov-pairs-folder', 
               type = str, 
               default = LOWCOV_PAIRS_SUBFOLDER,
               help = 'The name of the folder for lowfreq pairsam files.')
@click.option('--lowcov-sam-folder', 
               type = str, 
               default = LOWCOV_SAM_SUBFOLDER,
               help = 'The name of the folder for lowfreq sam files.')
@click.option('--highcov-pairs-folder', 
               type = str, 
               default = HIGHCOV_PAIRS_SUBFOLDER,
               help = 'The name of the folder for highcov pairsam files.')
@click.option('--highcov-sam-folder', 
               type = str, 
               default = HIGHCOV_SAM_SUBFOLDER,
               help = 'The name of the folder for highcov sam files.')
@click.option('--pairsam-ext', 
               default = PAIRS_EXT,
               type = str, 
               help = 'The file extension of the pairsam files.')
@click.option('--sam-ext', 
               default = SAM_EXT, 
               type = str, 
               help = 'The file extension of the sam files.')
@click.option('--stats-ext', 
               default = STATS_EXT, 
               type = str, 
               help = 'The file extension of the stats files.')
@click.option('--sort-threads', 
               default = 8, 
               type = int, 
               help = 'The amount of threads used by the sort command. Default: 8')
@click.option('--no-filterbycov', 
               is_flag = True, 
               help = 'Do not do the pairtools multifiler step after pairtools dedup. If set the --dups-pairs-folder option is ignored.')
def main(sam_folder, output_folder, output_sam_folder, chrom_sizes_file, assembly_name, working_dir, parse_stats_folder, dedup_stats_folder,
         filterbycov_stats_folder, unmapped_pairs_folder, unmapped_sam_folder, dups_pairs_folder, dups_sam_folder,
         lowcov_pairs_folder, lowcov_sam_folder, highcov_pairs_folder, highcov_sam_folder,
         pairsam_ext, sam_ext, stats_ext, sort_threads, no_filterbycov):
    '''Create the pairtools command for the snHiC pipeline.

    SAM_FOLDER: The directory with the sam files.

    OUTPUT_FOLDER: The folder with the resulting pairtools files.
    
    CHROM_SIZES_FILE: Chromosome order used to flip interchromosomal mates: path to a chromosomes file (e.g. UCSC chrom.sizes or similar) whose first column lists scaffol names. Any scaffolds not listed will be ordered lexicographically following the names provided.

    '''
    overwrite = False

    if working_dir is None:
        working_dir = os.getcwd()

    # the directories
    working_dir = os.path.abspath(working_dir)
    sam_folder = get_full_dir(sam_folder, working_dir)
    chrom_sizes_file = get_full_dir(chrom_sizes_file, working_dir)
    parse_stats_folder = get_full_dir(parse_stats_folder, working_dir)
    dedup_stats_folder = get_full_dir(dedup_stats_folder, working_dir)
    filterbycov_stats_folder = get_full_dir(filterbycov_stats_folder, working_dir)
    unmapped_pairs_folder = get_full_dir(unmapped_pairs_folder, working_dir)
    unmapped_sam_folder = get_full_dir(unmapped_sam_folder, working_dir)
    dups_pairs_folder = get_full_dir(dups_pairs_folder, working_dir)
    dups_sam_folder = get_full_dir(dups_sam_folder, working_dir)
    highcov_pairs_folder = get_full_dir(highcov_pairs_folder, working_dir)
    highcov_sam_folder = get_full_dir(highcov_sam_folder, working_dir)
    output_folder = get_full_dir(output_folder, working_dir)
    output_sam_folder = get_full_dir(output_sam_folder, working_dir)

    # the file extensions
    parse_stats_ext = PARSE_EXT.format(stats_ext)
    dedup_stats_ext = DEDUP_EXT.format(stats_ext)
    filterbycov_stats_ext = FILTERBYCOV_EXT.format(stats_ext)
    unmapped_pairs_ext = UNMAPPED_EXT.format(pairsam_ext)
    unmapped_sam_ext = UNMAPPED_EXT.format(sam_ext)
    dups_pairs_ext = DUPS_EXT.format(pairsam_ext)
    dups_sam_ext = DUPS_EXT.format(sam_ext)
    lowcov_pairs_ext = LOWCOV_EXT.format(pairsam_ext)
    lowcov_sam_ext = LOWCOV_EXT.format(sam_ext)
    highcov_pairs_ext = HIGHCOV_EXT.format(pairsam_ext)
    highcov_sam_ext = HIGHCOV_EXT.format(sam_ext)

    if no_filterbycov:
        output_ext = dups_pairs_ext
        output_sam_ext = dups_sam_ext
        
    else:
        output_ext = lowcov_pairs_ext
        output_sam_ext = lowcov_sam_ext

    cmds = list()
    sam_files = glob.glob(os.path.join(sam_folder, '*{}'.format(sam_ext)))
    sam_files.sort()
    
    for sam_file in sam_files:
        (_, sam_file_name) = os.path.split(sam_file)

        output_file = os.path.join(output_folder, sam_file_name.replace(sam_ext, output_ext))
        output_sam_file = os.path.join(output_sam_folder, sam_file_name.replace(sam_ext, output_sam_ext))

        parse_stats_file = os.path.join(parse_stats_folder, sam_file_name.replace(sam_ext, parse_stats_ext))
        dedup_stats_file = os.path.join(dedup_stats_folder, sam_file_name.replace(sam_ext, dedup_stats_ext))
        filterbycov_stats_file = os.path.join(filterbycov_stats_folder, sam_file_name.replace(sam_ext, filterbycov_stats_ext))
        unmapped_pairs_file = os.path.join(unmapped_pairs_folder, sam_file_name.replace(sam_ext, unmapped_pairs_ext))
        unmapped_sam_file = os.path.join(unmapped_sam_folder, sam_file_name.replace(sam_ext, unmapped_sam_ext))
        dups_pairs_file = os.path.join(dups_pairs_folder, sam_file_name.replace(sam_ext, dups_pairs_ext))
        dups_sam_file = os.path.join(dups_sam_folder, sam_file_name.replace(sam_ext, dups_sam_ext))
        
        highcov_pairs_file = os.path.join(highcov_pairs_folder, sam_file_name.replace(sam_ext, highcov_pairs_ext))
        highcov_sam_file = os.path.join(highcov_sam_folder, sam_file_name.replace(sam_ext, highcov_sam_ext))

        if not overwrite and os.path.exists(output_file):
            continue

        if no_filterbycov:
            cmd = CMD_DEDUP.format(parse_stats_file = parse_stats_file, assembly_name = assembly_name, chroms_path = chrom_sizes_file, sam_file = sam_file, threads = sort_threads, unmapped_pairs_file = unmapped_pairs_file, unmapped_sam_file = unmapped_sam_file, dedup_stats_file = dedup_stats_file, dups_pairs_file = dups_pairs_file, dups_sam_file = dups_sam_file, output_file = output_file, output_sam_file = output_sam_file)
        else:
            cmd = CMD_FILTERBYCOV.format(parse_stats_file = parse_stats_file, assembly_name = assembly_name, chroms_path = chrom_sizes_file, sam_file = sam_file, threads = sort_threads, unmapped_pairs_file = unmapped_pairs_file, unmapped_sam_file = unmapped_sam_file, dedup_stats_file = dedup_stats_file, dups_pairs_file = dups_pairs_file, dups_sam_file = dups_sam_file, filterbycov_stats_file = filterbycov_stats_file, highcov_pairs_file = highcov_pairs_file, highcov_sam_file = highcov_sam_file, output_file = output_file, output_sam_file = output_sam_file)

        cmds.append(cmd)

    print('\n'.join(cmds))
    sys.stderr.write('{} commands created\n'.format(len(cmds)))

def get_full_dir(dir, working_dir):
    if dir == '.':
        folder = os.path.abspath(dir)

    elif dir[:2] == './':
        folder = os.path.abspath(dir) 

    else:
        folder = os.path.join(working_dir, dir)

    if not os.path.exists(folder):
        os.makedirs(folder)

    return folder


if __name__ == "__main__":
    main()
