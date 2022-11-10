#! /usr/bin/env python

__version__ = '0.0.1-dev'

import sys
import glob
import click
import os.path

CMD = 'bwa mem -SP -t {threads} {reference_fasta} {reads1_fastq} {read2_fastq} > {output_sam}'

SAM_EXT = '.sam'
FASTQ_EXT = '.fastq'
FASTQ_PAIRED1_EXT = '.1{fastq_ext}'
FASTQ_PAIRED2_EXT = '.2{fastq_ext}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.argument('fastq_folder',
              type = str, 
              required = True)
@click.argument('reference_genome',
              type = str, 
              required = True)
@click.argument('sam-folder',
              type = str, 
              required = True)
@click.option('--fastq-ext', 
               default = FASTQ_EXT,
               type = str, 
               show_default = True,
               help = 'The file extension of the fastq files.')
@click.option('--sam-ext', 
               default = SAM_EXT, 
               type = str, 
               show_default = True,
               help = 'The file extension of the sam files.')
@click.option('--threads', 
               default = 20, 
               type = int, 
               show_default = True,
               help = 'The amount of threads used by bwa mem.')
def main(fastq_folder, reference_genome, sam_folder, fastq_ext, sam_ext, threads):
    '''Create the bwa commands.

    FASTQ_FOLDER: The directory with the fastq files.

    REFERENCE_GENOME: A fastq file with the reference chromosomes.
    
    SAM_FOLDER: The output directory for the sam files.

    '''
    if not os.path.exists(sam_folder):
        os.makedirs(sam_folder)


    cmds = list()
    fastq1_files = os.path.join(fastq_folder, '*{}'.format(FASTQ_PAIRED1_EXT.format(fastq_ext = fastq_ext)))
    fastq1_file_names = glob.glob(fastq1_files)
    fastq1_file_names.sort()
    
    for fastq1_file in fastq1_file_names:
        (_, fastq1_file_name) = os.path.split(fastq1_file)

        fastq2_file = fastq1_file.replace(FASTQ_PAIRED1_EXT.format(fastq_ext = fastq_ext), FASTQ_PAIRED2_EXT.format(fastq_ext = fastq_ext))
        sam_file = os.path.join(sam_folder, fastq1_file_name.replace(FASTQ_PAIRED1_EXT.format(fastq_ext = fastq_ext), sam_ext))

        if os.path.exists(sam_file):
            continue

        cmd = CMD.format(reference_fasta = reference_genome, reads1_fastq = fastq1_file, read2_fastq = fastq2_file, threads = threads, output_sam = sam_file)
        cmds.append(cmd)

    print('\n'.join(cmds))
    sys.stderr.write('{} commands created\n'.format(len(cmds)))

if __name__ == "__main__":
    main()
