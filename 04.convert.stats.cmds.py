#! /usr/bin/env python

__version__ = '0.0.1-dev'

import os
import sys
import glob
import click

CMD = '{snhic_path}/bin/convert.stats.to.json.py --parse-stats {parse_stats} --dedup-stats {dedup_stats} --filterbycov-stats {filterbycov_stats} --output {output_file}'

PARSE_EXT = '.parse.stats'
DEDUP_EXT = '.dedup.stats'
FILTERBYCOV_EXT = '.filterbycov.stats'
OUTPUT_EXT = '.stats.json'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.argument('parse_folder',
              type = str, 
              required = True)
@click.argument('dedup_folder',
              type = str, 
              required = True)
@click.argument('filterbycov_folder',
              type = str, 
              required = True)
@click.argument('output_folder',
              type = str, 
              required = True)
@click.option('--parse-ext', 
               default = PARSE_EXT, 
               type = str, 
               help = 'The file extension of the pairsamtools parse statistics files. Default: .parse.stats')
@click.option('--dedup-ext', 
               default = DEDUP_EXT, 
               type = str, 
               help = 'The file extension of the pairsamtools dedup statistics files. Default: .dedup.stats')
@click.option('--filterbycov-ext', 
               default = FILTERBYCOV_EXT, 
               type = str, 
               help = 'The file extension of the input pairsamtools filterbycov files. Default: .filterbycov.stats')
@click.option('--output-ext', 
               default = OUTPUT_EXT, 
               type = str, 
               help = 'The file extension of the output json files. Default: .stats.json')
@click.option('--snhic-path',  
               type = str,
               default = None,
               help = 'The path to the snHi-C pipeline resources. If $SNHIC_PATH is set this is not necessary.')
def main(parse_folder, dedup_folder, filterbycov_folder, output_folder, parse_ext, dedup_ext, filterbycov_ext, output_ext, snhic_path):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    parse_files = glob.glob(os.path.join(parse_folder, '*{}'.format(parse_ext)))
    
    if snhic_path == None:
        snhic_path = os.environ['SNHIC_PATH']

    cmds = list()
    for parse_file in parse_files:
        (_, parse_file_name) = os.path.split(parse_file)

        
        dedup_file_name = parse_file_name.replace(parse_ext, dedup_ext)
        dedup_file = os.path.join(dedup_folder, dedup_file_name)

        filterbycov_file_name = parse_file_name.replace(parse_ext, filterbycov_ext)
        filterbycov_file = os.path.join(filterbycov_folder, filterbycov_file_name)

        output_file_name = parse_file_name.replace(parse_ext, output_ext)
        output_file = os.path.join(output_folder, output_file_name)

        cmds.append(CMD.format(snhic_path = snhic_path, parse_stats = parse_file, dedup_stats = dedup_file, filterbycov_stats = filterbycov_file, output_file = output_file))

    print('\n'.join(cmds))
    sys.stderr.write('{} commands created\n'.format(len(cmds)))

if __name__ == '__main__':
    main()