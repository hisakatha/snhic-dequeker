#! /usr/bin/env python

__version__ = '0.0.1-dev'


import os
import sys
import glob
import click

import cooler
import pandas as pd
from snHiC.lib import coolerAnalysisShared as css

OUTPUT_EXT = '.{resolution}{output_ext}'
FILE_MERGE = True
TAG = 'tag1'

COOL_EXT = '.cool'
MCOOL_EXT = '.mcool'

FILE_NAME = '*{cellID}*.{resolution}{file_ext}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--samples-file', 
               required = True,
               type = click.Path(exists = True, readable = True), 
               help = 'The file with the sample meta data.')
@click.option('--cooler-dir', 
               required = True,
               type = click.Path(exists = True, readable = True),
               help = 'The directory with the .cool and .mcool files.')
@click.option('--output-dir', 
               required = True,
               type = click.Path(exists = True, readable = True),
               help = 'The directory for the merged .cool and .mcool files')
@click.option('--resolution', 
               default = 10000, 
               type = int, 
               help = 'The resolution of the .cool files that will be analzed.')
@click.option('--tags', 
               default = TAG,
               type = str, 
               help = 'Comma seperated lists of the tags used to group the samples. Default: {}'.format(TAG))
@click.option('--output-ext', 
               default = COOL_EXT,
               show_default = True,
               type = str, 
               help = 'The cooler file extension.')
def main(samples_file, cooler_dir, output_dir, resolution, tags, output_ext):
    '''Prepare the data for future use in generating the figures.
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

    grouped_dict = df.groupby(tag_names).groups
    #print(grouped_dict.keys())
    for key in grouped_dict.keys():
        if isinstance(key, str) or isinstance(key, int):
            output_cooler_file_name = '{}{}'.format(key, OUTPUT_EXT.format(resolution = resolution, output_ext = output_ext)) 
            #output_mcool_file_name = '{}{}'.format(key, OUTPUT_EXT.format(resolution = resolution, output_ext = mcool_ext))  
        else:
            arr = list()
            for k in key:
                if not k == '':
                    arr.append(k)
                    
            output_cooler_file_name = '{}{}'.format('_'.join([str(e) for e in arr]), OUTPUT_EXT.format(resolution = resolution, output_ext = output_ext))
            #output_mcool_file_name = '{}{}'.format('_'.join([str(e) for e in arr]), OUTPUT_EXT.format(resolution = resolution, output_ext = mcool_ext))
            
        output_cooler_file = os.path.join(output_dir, output_cooler_file_name)
        #print(output_cooler_file)
        #output_mcool_file = os.path.join(output_dir, output_mcool_file_name)
        if os.path.exists(output_cooler_file):
            continue
        
        merge_cooler_files = list()
        #merge_mcool_files = list()
        cellIDs = df.ix[grouped_dict[key]].cellID.values
        #print(cellIDs)
        for cellID in cellIDs:
            cooler_file_name = FILE_NAME.format(cellID = cellID, resolution = resolution, file_ext = output_ext)
            #mcool_file_name = FILE_NAME.format(cellID = cellID, resolution = resolution, file_ext = mcool_ext)

            cooler_file = os.path.join(cooler_dir, cooler_file_name)
            #mcool_file = os.path.join(cooler_dir, mcool_file_name)

            _cooler_files = glob.glob(cooler_file)
            #_mcool_files = glob.glob(mcool_file)

            if len(_cooler_files) > 0:
                merge_cooler_files.append(_cooler_files[0])

            #if len(_mcool_files) > 0:
            #    merge_mcool_files.append(_mcool_files[0])


        #print(output_mcool_file)
        #print(merge_mcool_files)
        #sys.exit()
        #print(merge_cooler_files)
        if not len(merge_cooler_files) == 0:
            css.mergeCoolerFiles(merge_cooler_files, output_cooler_file)

        print('{} files merged'.format(len(merge_cooler_files)))

        #css.mergeCoolerFiles(merge_mcool_files, output_mcool_file)
        #print('{} files merged'.format(len(merge_mcool_files)))

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

if __name__ == "__main__":
    main()
