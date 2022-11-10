#! /usr/bin/env python

__version__ = '0.0.1-dev'
import os
import sys
import glob
import click
import os.path

CMD_PE = 'samtools fastq -1 {output_file_1} -2 {output_file_2} --threads {threads} {input_file}'
CMD_SE = 'samtools fastq -0 {output_file} --threads {threads} {input_file}'
BAM_EXT = '.bam'
FASTQ_EXT = '.fastq'

FASTQ_PAIRED_EXT = '.{paired_end}{fastq_ext}'
BAM_FILE_NAMES = '*{bam_ext}'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.argument('bam_folder',
              type = str, 
              required = True)
@click.argument('fastq_folder',
              type = str, 
              required = True)
@click.option('--fastq-ext', 
               default = FASTQ_EXT,
               type = str, 
               help = 'The file extension of the fastq files.')
@click.option('--bam-ext', 
               default = BAM_EXT, 
               type = str, 
               help = 'The file extension of the bam files.')
@click.option('--threads', 
               default = 4, 
               type = int, 
               help = 'The amount of threads for samtools.')
@click.option('--pe/--se', 
               default = False, 
               show_default = True,
               help = 'Whether the sequencing was pair-ended oir single-ended.')
def main(bam_folder, fastq_folder, fastq_ext, bam_ext, threads, paired_ended):
    '''Prepare the data for future use in generating the figures.

    BAM_FOLDER: The directory with the bam files.

    FASTQ_FOLDER: The output directory for the fastq files.
    
    '''
    if not os.path.exists(fastq_folder):
        os.makedirs(fastq_folder)

    bam_files = os.path.join(bam_folder, BAM_FILE_NAMES.format(bam_ext = bam_ext))

    cmds = list()
    for bam_file in glob.glob(bam_files):
        (_, bam_file_name) = os.path.split(bam_file)
        fastq_file1 = os.path.join(fastq_folder, bam_file_name.replace(bam_ext, FASTQ_PAIRED_EXT.format(paired_end = 1, fastq_ext = fastq_ext)))
        fastq_file2 = os.path.join(fastq_folder, bam_file_name.replace(bam_ext, FASTQ_PAIRED_EXT.format(paired_end = 2, fastq_ext = fastq_ext)))

        if os.path.exists(fastq_file1) and os.path.exists(fastq_file2):
            continue

        cmd = CMD.format(output_file_1 = fastq_file1, output_file_2 = fastq_file2, threads = threads - 1, input_file = bam_file)
        cmds.append(cmd)

    print('\n'.join(cmds))
    sys.stderr.write('{} commands created\n'.format(len(cmds)))

if __name__ == "__main__":
    main()
