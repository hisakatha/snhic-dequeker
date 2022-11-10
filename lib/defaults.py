# system defaults
TAG = 'tag1'
VERSION = '0.0.1-dev'
RESOLUTION = 10000
SADDLE_DIM = 5
SMOOTH_SIGMA = 2
PLOTS_PER_PAGE = 4
INCLUDE_COLUMN = 'include'
SAVE_SOURCE_DATA = True
CONTEXT_SETTINGS = {'help_option_names': ['-h', '--help', '--rhonda'],}
MIN_EXP_LOOP_SIZE = 0
MAX_EXP_LOOP_SIZE = 1e6
METAKEYS = ['cis', 'cis_1kb+', 'cis_20kb+', 'trans', 'total', 'total_nodups', 'files_merged', 'avg_loop_size']


# file extension defaults
BAM_EXT = '.bam'
PDF_EXT = '.pdf'
SAM_EXT = '.sam'
FASTQ_EXT = '.fastq'
PAIRS_EXT = '.pairs'
STATS_EXT = '.stats'
COOLER_EXT = '.cool'
PARSE_EXT = '.parse{}'
DEDUP_EXT = '.dedup{}'
FILTERBYCOV_EXT = '.filterbycov{}'
UNMAPPED_EXT = '.unmapped{}'
DUPS_EXT = '.dups{}'
LOWCOV_EXT = '.lowcov{}'
HIGHCOV_EXT = '.highcov{}'
PARSE_STATS_EXT = '.parse.stats'
DEDUP_STATS_EXT = '.dedup.stats'
FILTERBYCOV_STATS_EXT = '.filterbycov.stats'
FINAL_STATS_EXT = '.stats.txt'
OUTPUT_STATS_EXT = '.stats.json'

# pickle file names
SADDLE_PKL_FILE = 'Saddle_{}.pkl'
PILEUP_LOOPS_PKL_FILE = 'PileUp_loops_{}.pkl'
SCALINGCOOLER_PKL_FILE = 'ScalingCooler_{}.pkl'
REGIONAVG_DOMAINS_PKL_FILE = 'RegionAvg_domains_{}.pkl'

# subfolders names defaults
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

