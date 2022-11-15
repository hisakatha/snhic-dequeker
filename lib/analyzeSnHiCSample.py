import os
import cooler
import pickle
import bioframe
import numpy as np
import pandas as pd
from scipy.linalg import toeplitz
from sinkhorn_knopp import sinkhorn_knopp as skp
from lib import coolerAnalysisShared as css

from cooltools.api import expected
from cooltools.api import eigdecomp

OVERWRITE = False

SADDLE_FILE_NAME = 'Saddle_{}.pkl'
PILE_UP_FILE_NAME = 'PileUp_{}_{}.pkl'
REGION_AVG_FILE_NAME = 'RegionAvg_{}_{}.pkl'
SCALING_COOL_FILE_NAME = 'ScalingCooler_{}.pkl'

# make plots for figure 1
def deploySample(data_file, save_folder, genome_path, scaling_ref_file, loopAndDomainDict, saddle_plot, overwrite = OVERWRITE):
    if not os.path.exists(save_folder):
        print('creating: {}'.format(save_folder))
        os.makedirs(save_folder)

    assert(os.path.exists(save_folder))

    if scaling_ref_file == None:
        print('warning: reference cooler file is not available')
    elif not os.path.exists(scaling_ref_file):
        print('warning: reference cooler file {} does not exist'.format(scaling_ref_file))

    ####################
    # needs to be removed
    #saddle_plot = True

    print("processing {}".format(data_file))
    doOnePileUp(data_file, save_folder, loopAndDomainDict, overwrite)
    if saddle_plot:
        doOneSaddle(data_file, save_folder, genome_path, overwrite)

    doOneRescaledPileUp(data_file, save_folder, loopAndDomainDict, overwrite)
    doOneCoolerScaling(data_file, save_folder, scaling_ref_file, overwrite)

def deploySamples(existing_csv, save_folder, scaling_ref_file, loopAndDomainDict, overwrite = OVERWRITE):
    good_df = pd.read_csv(existing_csv)
    lof = list(zip(good_df.filepath.values, good_df.filenames.values))

    for (datapath, filename) in lof:
        data_file = os.path.join(datapath, filename)
        deploySample(data_file, save_folder, scaling_ref_file, loopAndDomainDict, overwrite)

########### Aggregate Analyses ###############
def doOnePileUp(data_file, save_folder, loopAndDomainDict, overwrite):
    
    loopTypes = ['loops']
    separations = [100000, 100000000]
    separations = [100000, 150000, 250000, 500000]
    padsize = 10 # 120 kb if resolution = 10000
    
    stored_maps = dict()     
    for loopKeys in loopTypes:
        
        (_, file_name) = os.path.split(data_file)
        fname = os.path.join(save_folder, PILE_UP_FILE_NAME.format(loopKeys, file_name))
        
        # don't re-do if overwrite = false
        if os.path.exists(fname) and not overwrite:
            print('warning: {} exists and cannot overwrite\nskipping'.format(fname))
            return 
        
        loopPositions = loopAndDomainDict[loopKeys]
        loopPositions['sep'] = loopPositions['ends'] - loopPositions['starts'] 
        
        for (_, lo, hi) in zip(range(len(separations[0:-1])), separations[0:-1], separations[1:]):

            goodLoops = (loopPositions['sep'] >= lo) & (loopPositions['sep'] < hi)
            these_loops = loopPositions.loc[goodLoops]

            (M, C, m_count, c_count) = css.averageLoopsWithControl(these_loops, data_file, numShifted = 100, pad = padsize)

            stored_maps['{0}-{1}bp'.format(lo, hi)] = (file_name, loopKeys, M, C, m_count, c_count)
        
        print('writing {}'.format(fname))
        pickle.dump(stored_maps, open(fname, 'wb'))
        print('done writing {}'.format(fname))
    
def doOneRescaledPileUp(data_file, save_folder, loopAndDomainDict, overwrite): 
    #print(data_file)
    loopTypes = ['domains']
    for loopKeys in loopTypes:
        (_, file_name) = os.path.split(data_file)
        fname = os.path.join(save_folder, REGION_AVG_FILE_NAME.format(loopKeys, file_name))

        # don't re-do if overwrite = false
        if os.path.exists(fname) and not overwrite:
            print('warning: {} exists and cannot overwrite\nskipping'.format(fname))
            return  
        
        domainPositions = loopAndDomainDict[loopKeys]
        maps = css.averageDomains(data_file, minSize = 100000, maxSize = 1000000, domainPositions = domainPositions)
        
        print('writing {}'.format(fname))
        pickle.dump(maps, open(fname, 'wb'))
        print('done writing {}'.format(fname))

def doOneCoolerScaling(data_file, save_folder, scaling_ref_file, overwrite):
    (_, file_name) = os.path.split(data_file)
    fname = os.path.join(save_folder, SCALING_COOL_FILE_NAME.format(file_name))
    
    # don't re-do if overwrite = false
    if os.path.exists(fname) and not overwrite:
        print('warning: {} exists and cannot overwrite\nskipping'.format(fname))
        return        
                         
    scaling_cooler = css.getCoolerScaling(data_file, reference_cooler_file = scaling_ref_file) #logFactor = 1.15   
    
    print('writing {}'.format(fname))
    pickle.dump(scaling_cooler, open(fname, 'wb'))
    print('done writing {}'.format(fname))

def doOneSaddle(data_file, save_folder, genome_path, overwrite):
    (_, file_name) = os.path.split(data_file)
    fname = os.path.join(save_folder, SADDLE_FILE_NAME.format(file_name))
    
    # don't re-do if overwrite = false
    if os.path.exists(fname) and not overwrite:
        print('warning: {} exists and cannot overwrite\nskipping'.format(fname))
        return   

    saddles = css.doCoolerSaddle(data_file, genome_path)

    print('writing {}'.format(fname))
    pickle.dump(saddles, open(fname, 'wb'))
    print('done writing {}'.format(fname))     


def compute_eigenvectors(genome_chromsizes, cool_file, long_name, bins, resolution, data_folder, balance = False, overwrite = OVERWRITE):
    #######################################################################
    # compute the eigenvectors
    
    lam_csv_file_name = '{name}.{res}.eigs.cis.lam.txt'.format(name = long_name, res = resolution)
    lam_csv_file = os.path.join(data_folder, lam_csv_file_name)
    eigs_csv_file_name = '{name}.{res}.eigs.cis.vecs.txt'.format(name = long_name, res = resolution)
    eigs_csv_file = os.path.join(data_folder, eigs_csv_file_name)
    
    if overwrite or not (os.path.exists(lam_csv_file) and os.path.exists(eigs_csv_file)):
        if balance:
            (lam, eigs) = eigdecomp.cooler_cis_eig(cool_file, bins, n_eigs = 3, phasing_track_col = 'GC', sort_metric='var_explained')
        else:
            (lam, eigs) = eigdecomp.cooler_cis_eig(cool_file, bins, n_eigs = 3, phasing_track_col = 'GC', balance = False, ignore_diags = 2, sort_metric='var_explained')
        
        lam.to_csv(lam_csv_file, sep='\t')
        eigs.to_csv(eigs_csv_file, sep='\t', index = False)
    else:
        print('warning: {} exists and was not overwritten'.format(os.path.split(lam_csv_file)[1]))
        print('warning: {} exists and was not overwritten'.format(os.path.split(eigs_csv_file)[1]))
    
    # Save bigwig track
    bigwig_track_file_name = '{name}.{res}.eigs.cis.vecs.E1.bw'.format(name = long_name, res = resolution)
    bigwig_track_file = os.path.join(data_folder, bigwig_track_file_name)
    if overwrite or not os.path.exists(bigwig_track_file):
        bioframe.to_bigwig(eigs, genome_chromsizes, bigwig_track_file, 'E1')
    else:
        print('warning: {} exists and was not overwritten'.format(os.path.split(bigwig_track_file)[1]))


    
def get_gc_content(genome_chromsizes, resolution, genome_fasta):
    #######################################################################
    # get the genome GC content
    #print('binning')
    bins = cooler.binnify(genome_chromsizes, resolution)
    #print('loading')
    fasta_records = bioframe.load_fasta(genome_fasta)
    #print('gc content')
    bins['GC'] = bioframe.tools.frac_gc(bins, fasta_records)
    #print(bins.head())

    return bins


def generate_cis_trans_stats(genome_chromsizes, cool_file, long_name, resolution, data_folder, balance = False, overwrite = OVERWRITE):
    ####################################################################
    # generate the expected cis and trans stats
    chroms = genome_chromsizes.index
    supports = [(chrom, 0, genome_chromsizes[chrom]) for chrom in chroms]
    
    expected_cis_file_name = '{name}.{res}.expected.cis.tsv'.format(name = long_name, res = resolution)
    expected_cis_file = os.path.join(data_folder, expected_cis_file_name)
    #expected_trans_file_name = '{name}.{res}.expected.trans.tsv'.format(name = long_name, res = resolution)
    #expected_trans_file = os.path.join(data_folder, expected_trans_file_name)

    
    transforms = dict()
    if overwrite or not os.path.exists(expected_cis_file):
        if balance:
            transforms['balance'] = lambda p: p['count'] * p['weight1'] * p['weight2']
            tables = expected.diagsum(cool_file, supports, transforms = transforms, chunksize = 10000000, ignore_diags = 2)
        else:
            tables = expected.diagsum(cool_file, supports, transforms = transforms, chunksize = 10000000, ignore_diags = 2)
        
        cis_exp = pd.concat([tables[support] for support in supports], keys = [support[0] for support in supports], names = ['chrom'])
        
        if balance:
            cis_exp['balance.avg'] = cis_exp['balance.sum'] / cis_exp['n_valid']
        
        cis_exp['count.avg'] = cis_exp['count.sum'] / cis_exp['n_valid']
        print(cis_exp)
        cis_exp.to_csv(expected_cis_file, sep='\t')
    else:
        print('warning: {} exists and was not overwritten'.format(os.path.split(expected_cis_file)[1]))
    
    """
    if overwrite or not os.path.exists(expected_trans_file):
        #records = expected.blocksum_pairwise(cool_file, supports, transforms = transforms, chunksize = 10000000)
        records = expected.blocksum_pairwise(cool_file, supports, transforms = dict(), chunksize = 10000000)
        trs_exp = pd.DataFrame([{'chrom1': s1[0], 'chrom2': s2[0], **rec} for (s1, s2), rec in records.items()], columns = ['chrom1', 'chrom2', 'n_valid', 'count.sum', 'balanced.sum'])
        trs_exp.to_csv(expected_trans_file, sep = '\t')
    else:
        print('warning: {} exists and was not overwritten'.format(os.path.split(expected_trans_file)[1]))
    """

def get_chrom_sizes(chrom_sizes_file):
    f = open(chrom_sizes_file)
    d = dict()
    for line in f:
        line = line.strip()
        try:
            (chrom_name, chrom_size) = line.split('\t')
            chrom_size = int(chrom_size)
            d[chrom_name] = chrom_size
        except:
            pass

    chrom_sizes = pd.Series(d, d.keys())
    chrom_sizes.name = 'length'

    return chrom_sizes

def make_cis_obsexp_fetcher(clr, expected):
    """
    Construct a function that returns intra-chromosomal OBS/EXP.

    Parameters
    ----------
    clr : cooler.Cooler
        Observed matrix.
    expected : DataFrame
        Diagonal summary statistics for each chromosome.
    name : str
        Name of data column in ``expected`` to use.

    Returns
    -------
    getexpected(chrom, _). 2nd arg is ignored.

    """
    expected, name = expected
    expected = {k: x.values for k, x in expected.groupby('chrom')[name]}
    return lambda chrom, _: (clr.matrix(balance = False).fetch(chrom) / toeplitz(expected[chrom]))

def balance_saddle_data(saddle_data):
    matrix_size = len(saddle_data)
    A = list()
    for i in range(1, matrix_size - 1):
        A.append(saddle_data[i][1:matrix_size - 1])

    sk = skp.SinkhornKnopp()
    P_ds = sk.fit(A)

    B = np.zeros((matrix_size, matrix_size), dtype = float)
    B.fill(np.nan)

    """
    for (i, saddle_arr) in enumerate(saddle_data):
        for (j, element) in enumerate(saddle_arr):
           B[i][j] = element

    print(B)
    """
    
    for (i, P_ds_arr) in enumerate(P_ds):
        for (j, element) in enumerate(P_ds_arr):
            B[i + 1][j + 1] = element


    return B


    





    
