#! /usr/bin/env python

__version__ = '0.0.1-dev'
import glob
import click
import os.path
import pandas as pd

DATA_SUB_DIR = 'data'
ANALYSES_SUB_DIR = 'analyses'

DATA_PDF_SUB_DIR = 'data.pdf'
DATA_FRAMES_SUB_DIR = 'data.frames'
DATA_PROCESSED_SUB_DIR = 'data.processed'

SCRIPTS_SUB_DIR = 'scripts'
SLURM_SCRIPTS_SUB_DIR = 'slurm.scripts'
COOLERS_MERGED_SUB_DIR = 'coolers.merged'
TAD_INSUL_SUB_DIR = 'tad.insul'
CLUSTER_OUTPUT_SUB_DIR = 'log'

BAM_SUB_DIR = 'bam'
SAM_SUB_DIR = 'sam'
FASTQ_SUB_DIR = 'fastq'
PAIRS_LOWCOV_SUB_DIR = 'pairs.lowcov'
SAM_LOWCOV_SUB_DIR = 'sam.lowcov'
STATS_PARSE_SUB_DIR = 'stats.parse'
STATS_DEDUP_SUB_DIR = 'stats.dedup'
STATS_FILTERBYCOV_SUB_DIR = 'stats.filterbycov'
STATS_SUB_DIR = 'stats'
COOLER_SUB_DIR = 'cooler'
BOOTSTRAP_SUB_DIR = 'bootstraps'
FINAL_STATS_SUB_DIR = 'stats'

ZOOMIFY_INPUT_EXT = '.cool'
TMPL_FILE_NAMES = '*.tmpl'
SCRIPTS_FILE_NAME = '{prefix}.sh'

CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
@click.version_option(version = __version__)
@click.command(context_settings = CONTEXT_SETTINGS)
@click.option('--xlsx-file', 
               type = str, 
               required = True,
               help = 'The xlsx file with multiple project setups.')
@click.option('--scripts-tmpl-folder', 
               type = str, 
               required = True,
               help = 'The folder with the scripts templates.')
@click.option('--slurm-tmpl-folder', 
               type = str, 
               required = True,
               help = 'The folder with the slurm templates.')
def main(xlsx_file, scripts_tmpl_folder, slurm_tmpl_folder):
    '''
    
    '''

    df = pd.read_excel(xlsx_file, sep = ',')
    df = df[df['include'] == True]

    project_id = df.projectID.values
    analysis_name = df.analysisName.values
    paired_conditions = df.pairedConditions.values
    project_folder = df.projectFolder.values
    samples_file = df.samplesFile.values
    loops_file = df.loopsFile.values

    domains_file = df.domainsFile.values
    scaling_ref_file = df.scalingRefFile.values
    max_exp_loop_size = df.maxLoopSize.values
    group_tags = df.groupTags.values
    aggr_resolution = df.aggrResolution.values
    cooler_resolutions = df.coolerResolutions.values
    zoomify_resolutions = df.zoomifyResolutions.values

    assembly_name = df.assemblyName.values
    genome_file = df.genomeFile.values
    genome_folder = df.genomeFolder.values
    chrom_sizes_file = df.chromSizesFile.values
    bwa_threads = df.bwaThreads.values
    sort_threads = df.sortThreads.values
    samtools_threads = df.samtoolsThreads.values
    tad_insul_diag = df.tadInsulationDiagonals.values
    tad_insul_diamond = df.tadInsulationDiamondSize.values
    tad_insul_dist = df.tadInsulationDistance.values
    
    plots_per_page = df.plotsPerPage.values
    saddle_plots = df.saddlePlots.values
    balance = df.balance.values
    by_chrom = df.byChrom.values

    bootstrap_threads = df.bootstrapThreads.values
    bs_loop_strength_iterations = df.bsLoopStrengthIterations.values
    bs_loop_size_iterations = df.bsLoopSizeIterations.values
    rand_loop_strength_iterations = df.randLoopStrengthIterations.values

    first_script = df.firstScript.values

    projects_settings = zip(project_id, analysis_name, paired_conditions, project_folder, samples_file, loops_file, domains_file, scaling_ref_file, group_tags, aggr_resolution, cooler_resolutions, zoomify_resolutions, assembly_name, genome_file, genome_folder, chrom_sizes_file, bwa_threads, sort_threads, samtools_threads, plots_per_page, saddle_plots, balance, by_chrom, first_script, max_exp_loop_size, tad_insul_diag, tad_insul_diamond, tad_insul_dist, bootstrap_threads, bs_loop_strength_iterations, bs_loop_size_iterations, rand_loop_strength_iterations)
    
    #for project_settings in projects_settings:
    for project_setting in projects_settings:
        #(project_id, analysis_name, project_folder, samples_file, loops_file, domains_file, scaling_ref_file, group_tags, cooler_resolutions, zoomify_resolutions, assembly_name, genome_file, genome_folder, chrom_sizes_file, bwa_threads, sort_threads, samtools_threads, plots_per_page, saddle_plots, balance, by_chrom, first_script, max_exp_loop_size, tad_insul_diag, tad_insul_diamond, tad_insul_dist, bootstrap_threads, bs_loop_strength_iterations, bs_loop_size_iterations, rand_loop_strength_iterations) = project_settings
        (project_id, analysis_name, paired_conditions, project_folder, samples_file, loops_file, domains_file, scaling_ref_file, group_tags, aggr_resolution, cooler_resolutions, zoomify_resolutions, assembly_name, genome_file, genome_folder, chrom_sizes_file, bwa_threads, sort_threads, samtools_threads, plots_per_page, saddle_plots, balance, by_chrom, first_script, max_exp_loop_size, tad_insul_diag, tad_insul_diamond, tad_insul_dist, bootstrap_threads, bs_loop_strength_iterations, bs_loop_size_iterations, rand_loop_strength_iterations) = project_setting
        print('creating: {}'.format(analysis_name))
        (scripts_dir, slurm_dir, tmpl_vars) = create_project(project_id, analysis_name, paired_conditions, project_folder, samples_file, loops_file, domains_file, scaling_ref_file, group_tags, aggr_resolution, cooler_resolutions, zoomify_resolutions, assembly_name, genome_file, genome_folder, chrom_sizes_file, bwa_threads, sort_threads, samtools_threads, plots_per_page, saddle_plots, balance, by_chrom, max_exp_loop_size, tad_insul_diag, tad_insul_diamond, tad_insul_dist, bootstrap_threads, bs_loop_strength_iterations, bs_loop_size_iterations, rand_loop_strength_iterations)
        
        write_tmpl_files(scripts_tmpl_folder, scripts_dir, tmpl_vars, first_script)
        write_tmpl_files(slurm_tmpl_folder, slurm_dir, tmpl_vars, first_script)

def create_project(project_id, analysis_name, paired_conditions, project_folder, samples_file, loops_file, domains_file, scaling_ref_file, group_tags, aggr_resolution, cooler_resolutions, zoomify_resolutions, assembly_name, genome_file, genome_folder, chrom_sizes_file, bwa_threads, sort_threads, samtools_threads, plots_per_page, saddle_plots, balance, by_chrom, max_exp_loop_size, tad_insul_diag, tad_insul_diamond, tad_insul_dist, bootstrap_threads, bs_loop_strength_iterations, bs_loop_size_iterations, rand_loop_strength_iterations):
    if not os.path.exists(project_folder):
        os.makedirs(project_folder)

    data_dir = _makedirs(os.path.join(project_folder, DATA_SUB_DIR))
    analyses_sub_dir = _makedirs(os.path.join(project_folder, ANALYSES_SUB_DIR))
    analysis_dir = _makedirs(os.path.join(analyses_sub_dir, analysis_name))

    data_pdf_dir = _makedirs(os.path.join(analysis_dir, DATA_PDF_SUB_DIR))
    data_processed_dir = _makedirs(os.path.join(analysis_dir, DATA_PROCESSED_SUB_DIR))
    #data_frames_dir = _makedirs(os.path.join(analysis_dir, DATA_FRAMES_SUB_DIR))
    data_frames_dir = os.path.join(analysis_dir, DATA_FRAMES_SUB_DIR)
    scripts_dir = _makedirs(os.path.join(analysis_dir, SCRIPTS_SUB_DIR))
    slurm_scripts_dir = _makedirs(os.path.join(analysis_dir, SLURM_SCRIPTS_SUB_DIR))
    coolers_merged_dir = _makedirs(os.path.join(analysis_dir, COOLERS_MERGED_SUB_DIR))
    tad_insul_dir = _makedirs(os.path.join(analysis_dir, TAD_INSUL_SUB_DIR))
    cluster_output_dir = _makedirs(os.path.join(analysis_dir, CLUSTER_OUTPUT_SUB_DIR))
    cooler_dir = os.path.join(analysis_dir, COOLER_SUB_DIR)
    bootstrap_dir = os.path.join(analysis_dir, BOOTSTRAP_SUB_DIR)
    final_stats_dir = _makedirs(os.path.join(analysis_dir, FINAL_STATS_SUB_DIR))
    
    bam_dir = _makedirs(os.path.join(data_dir, BAM_SUB_DIR))
    sam_dir = os.path.join(data_dir, SAM_SUB_DIR)
    fastq_dir = os.path.join(data_dir, FASTQ_SUB_DIR)
    pairs_lowcov_dir = os.path.join(data_dir, PAIRS_LOWCOV_SUB_DIR)
    sam_lowcov_dir = os.path.join(data_dir, SAM_LOWCOV_SUB_DIR)
    stats_parse_dir = os.path.join(data_dir, STATS_PARSE_SUB_DIR)
    stats_dedup_dir = os.path.join(data_dir, STATS_DEDUP_SUB_DIR)
    stats_filterbycov_dir = os.path.join(data_dir, STATS_FILTERBYCOV_SUB_DIR)
    stats_dir = os.path.join(data_dir, STATS_SUB_DIR)
    
    
    #bwa_threads = bwa_threads - 1
    bootstrap_threads = bootstrap_threads
    resolution = cooler_resolutions.split(',')[0]

    tmpl_vars = dict()
    tmpl_vars['?BAM_DIR?'] = bam_dir
    tmpl_vars['?FASTQ_DIR?'] = fastq_dir
    tmpl_vars['?SLURM_DIR?'] = slurm_scripts_dir
    tmpl_vars['?BWA_THREADS?'] = int(bwa_threads - 1)
    tmpl_vars['?BWA_SLURM_THREADS?'] = int(bwa_threads)
    tmpl_vars['?SORT_THREADS?'] = int(sort_threads)
    tmpl_vars['?SAMTOOLS_THREADS?'] = int(samtools_threads)
    tmpl_vars['?PAIRED_CONDITIONS?'] = paired_conditions
    tmpl_vars['?SAMPLES_FILE?'] = samples_file
    tmpl_vars['?GENOME_FILE?'] = genome_file
    tmpl_vars['?SAM_DIR?'] = sam_dir
    tmpl_vars['?ASSEMBLY_NAME?'] = assembly_name
    tmpl_vars['?DATA_DIR?'] = data_dir
    tmpl_vars['?PAIRTOOLS_OUTPUT_PAIRS_DIR?'] = pairs_lowcov_dir
    tmpl_vars['?PAIRTOOLS_OUTPUT_SAM_DIR?'] = sam_lowcov_dir
    tmpl_vars['?STATS_PARSE_DIR?'] = stats_parse_dir
    tmpl_vars['?STATS_DEDUP_DIR?'] = stats_dedup_dir
    tmpl_vars['?STATS_FILTERBYCOV_DIR?'] = stats_filterbycov_dir
    tmpl_vars['?STATS_DIR?'] = stats_dir
    tmpl_vars['?AGGR_RESOLUTION?'] = aggr_resolution
    tmpl_vars['?CLOAD_RESOLUTIONS?'] = cooler_resolutions
    tmpl_vars['?COOLER_DIR?'] = cooler_dir
    tmpl_vars['?ZOOMIFY_INPUT_EXT?'] = '.{resolution}{zoomify_ext}'.format(resolution = resolution, zoomify_ext = ZOOMIFY_INPUT_EXT)
    tmpl_vars['?ZOOMIFY_RESOLUTIONS?'] = zoomify_resolutions
    tmpl_vars['?DATA_FRAMES_DIR?'] = data_frames_dir
    tmpl_vars['?RESOLUTION?'] = int(resolution)
    tmpl_vars['?TAGS?'] = group_tags
    tmpl_vars['?COOLERS_MERGED_DIR?'] = coolers_merged_dir
    tmpl_vars['?TAD_INSUL_DIR?'] = tad_insul_dir
    tmpl_vars['?BALANCE_OPTION?'] = '--balance' if balance else '--no-balance'
    tmpl_vars['?SADDLE_PLOTS_OPTION?'] = '--saddle-plots' if saddle_plots else '--no-saddle-plots'
    tmpl_vars['?BY_CHROM_OPTION?'] = '--by-chrom' if by_chrom else '--by-genome' 
    tmpl_vars['?GENOME_DIR?'] = genome_folder
    tmpl_vars['?LOOPS_FILE?'] = loops_file
    tmpl_vars['?DOMAINS_FILE?'] = domains_file
    tmpl_vars['?SCALING_REF_FILE?'] = scaling_ref_file
    tmpl_vars['?DATA_PROCESSED_DIR?'] = data_processed_dir
    tmpl_vars['?CHROM_SIZES_FILE?'] = chrom_sizes_file
    tmpl_vars['?DATA_PDFS_DIR?'] = data_pdf_dir
    tmpl_vars['?PLOTS_PER_PAGE?'] = int(plots_per_page)
    tmpl_vars['?FIGURE_NAME?'] = '{analysis_name}.{resolution}'.format(analysis_name = analysis_name, resolution = _convert_resolution(int(resolution)), balanced = 'balance' if balance else 'no-balance')
    tmpl_vars['?CLUSTER_OUTPUT_DIR?'] = cluster_output_dir
    tmpl_vars['?MAX_EXP_LOOP_SIZE?'] = int(max_exp_loop_size)
    tmpl_vars['?TAD_INSUL_DIAGONALS?'] = int(tad_insul_diag)
    tmpl_vars['?TAD_INSUL_DIAMOND_SIZE?'] = int(tad_insul_diamond)
    tmpl_vars['?TAD_INSUL_DISTANCE?'] = int(tad_insul_dist)
    tmpl_vars['?BOOTSTRAP_THREADS?'] = int(bootstrap_threads)
    tmpl_vars['?BOOTSTRAP_SLURM_THREADS?'] = int(bootstrap_threads) + 1
    tmpl_vars['?BS_LOOP_SIZE_ITERATIONS?'] = int(bs_loop_size_iterations)
    tmpl_vars['?BS_LOOP_STRENGTH_ITERATIONS?'] = int(bs_loop_strength_iterations)
    tmpl_vars['?RAND_LOOP_STRENGTH_ITERATIONS?'] = int(rand_loop_strength_iterations)
    tmpl_vars['?BOOTSTRAP_OUTPUT_DIR?'] = bootstrap_dir
    tmpl_vars['?FINAL_STATS_DIR?'] = final_stats_dir

    return (scripts_dir, slurm_scripts_dir, tmpl_vars)

def write_tmpl_files(tmpl_dir, scripts_dir, tmpl_vars, first_script):
    tmpl_files = glob.glob(os.path.join(tmpl_dir, TMPL_FILE_NAMES))
    tmpl_files.sort()
    for tmpl_file in tmpl_files:
        (_, tmpl_file_name) = os.path.split(tmpl_file)
        (tmpl_file_name_prefix, _) = os.path.splitext(tmpl_file_name)

        try:
            tmpl_file_pos = int(tmpl_file_name.split('.')[0])
            if not tmpl_file_pos >= first_script:
                continue
        except:
            pass

        #print(tmpl_file_name)
        
        scripts_file_name = SCRIPTS_FILE_NAME.format(prefix = tmpl_file_name_prefix)
        scripts_file = os.path.join(scripts_dir, scripts_file_name)

        f_tmpl = open(tmpl_file, 'r')
        tmpl = f_tmpl.read()

        for tmpl_var in tmpl_vars:
            tmpl = tmpl.replace(tmpl_var, str(tmpl_vars[tmpl_var]))
        
        f = open(scripts_file, 'w')
        f.write(tmpl)
        f.close()


def _makedirs(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

    return folder

def _convert_resolution(resolution):
    K = 1000
    M = 1000000
    G = 1000000000

    if resolution < K:
        return resolution
    elif K <= resolution < M:
        return "{:.0f}k".format(resolution / K)
    elif M <= resolution < G:
        return "{:.0f}M".format(resolution / M)
    else:
        return "{:.0f}G".format(resolution / G)

if __name__ == "__main__":
    main()