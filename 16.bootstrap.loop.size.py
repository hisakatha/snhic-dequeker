#! /usr/bin/env python
# (c) 2019. All Rights Reserved.
#       Sean Powell (sean.powell@imba.oeaw.ac.at)

__version__ = '0.0.1-dev'

import os
import sys
import click
import pandas as pd

OUTPUT_FILE_NAME = '{prefix}.loop.size.bs.txt'
COOLER_EXT = '.{resolution}{cooler_ext}'
TAG = 'tag1'
COOL_EXT = '.cool'

TMP_OPTION = '--tmp-dir {}'

CMD = '{snhic_path}/bin/bootstrap.loop.size.py --cooler-file-prefix {cooler_file_prefix} --sample-ids {sample_ids} {tmp_option} --cooler-dir {cooler_dir} --scaling-ref-file {scaling_ref_file} --output-file {output_file} --threads {threads} --iterations {iterations}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--samples-file', 
               type = click.Path(exists = True, resolve_path = True),
               required = True)
@click.option('--output-dir',
               type = click.Path(exists = True, resolve_path = True), 
               required = True)
@click.option('--scaling-ref-file', 
               required = True,
               type = click.Path(exists = True, resolve_path = True),
               help = 'A merged cool file (.mcool) with the scaling that can be used as a reference.')
@click.option('--tmp-dir',
              type = str, 
              help = 'The temporary directory used for the bootstrapping cooler files. If no directory is set the system tmp directory will bve used.')             
@click.option('--tags', 
               required = True,
               type = str, 
               help = 'The tag names seperated by a comma for the grouping of the samples.')
@click.option('--cooler-dir', 
               required = True,
               type = str, 
               help = 'The directory with the sample cooler files.')
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
@click.option('--snhic-path',  
               type = str,
               default = None,
               help = 'The path to the snHi-C pipeline resources. If $SNHIC_PATH is set this is not necessary.')
def main(samples_file, output_dir, scaling_ref_file, tmp_dir, tags, cooler_dir, cooler_ext, resolution, threads, iterations, snhic_path):
    '''Create bootstrapping command for the aggregate loop and domain analysis.

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

    if not tmp_dir == None:
        tmp_option = TMP_OPTION.format(tmp_dir)
    else:
        tmp_option = ''

    if snhic_path == None:
        snhic_path = os.environ['SNHIC_PATH']

    grouped_dict = df.groupby(tag_names).groups
    cmds = list()
    # create the prefix
    for key in grouped_dict.keys():
        if isinstance(key, str) or isinstance(key, int):
            cooler_file_prefix = key
            #cooler_file_name = '{}{}'.format(key, COOLER_EXT.format(resolution = resolution, cooler_ext = cooler_ext)) 
            #output_file_name = OUTPUT_FILE_NAME.format(prefix = key)
        else:
            arr = list()
            for k in key:
                if not k == '':
                    arr.append(k)

            cooler_file_prefix = '_'.join([str(e) for e in arr])

            #cooler_file_name = '{}{}'.format('_'.join([str(e) for e in arr]), COOLER_EXT.format(resolution = resolution, output_ext = cooler_ext))
            #output_file_name = OUTPUT_FILE_NAME.format(prefix = '_'.join([str(e) for e in arr]))
            
        #cooler_file = os.path.join(cooler_dir, cooler_file_name)
        output_file = os.path.join(output_dir, OUTPUT_FILE_NAME.format(prefix = cooler_file_prefix))

        cellIDs = df.ix[grouped_dict[key]].cellID.values
        cellIDs.sort()
        sample_ids = ','.join([str(cellID) for cellID in cellIDs])

        cmd = CMD.format(snhic_path = snhic_path, cooler_file_prefix = cooler_file_prefix, sample_ids = sample_ids, tmp_option = tmp_option, cooler_dir = cooler_dir, scaling_ref_file = scaling_ref_file, output_file = output_file, resolution = resolution, threads = threads, iterations = iterations)
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