# (c) 2019. All Rights Reserved.
# Code written by: Hugo Brandao (hbrandao@g.harvard.edu)
# Code written by: Maksim Imakaev (imakaev@mit.edu)
# Code written by: Sean Powell (sean.powell@imba.oeaw.ac.at)

import re
import six
import os
import sys
#import joblib
import pandas as pd
import numpy as np

from mirnylib.numutils import zoomArray, coarsegrain, ultracorrect, observedOverExpected, logbinsnew
from hiclib.hicShared import getResolution

import cooler
from functools import reduce
import pickle

from mirnylib.genome import Genome

# mem = joblib.Memory(".")  # cache


def getLoopsAndDomains(genome_path, loops_file, domains_file):
    mygen = Genome(genome_path, readChrms = ["#", "X"])

    loops_df = pd.read_csv(loops_file, sep="\t")
    loops = getLoops(mygen, loops_df)
    
    domains_df = pd.read_csv(domains_file, sep="\t")
    domains = getDomains(mygen, domains_df)
    
    CUTOFF = 80000
    domStartID = domains["chrms"].values * 500000000 + domains["starts"].values
    domEndID = domains["chrms"].values * 500000000 + domains["ends"].values
    loopStartID = loops["chrms"].values * 500000000 + loops["starts"].values
    loopEndID = loops["chrms"].values * 500000000 + loops["ends"].values

    mask = list()
    for st, end in zip(loopStartID, loopEndID):
        if (np.abs(domStartID - st) + np.abs(domEndID - end)).min() < CUTOFF:
            mask.append(True)
        else:
            mask.append(False)

    domainsLoops = loops.ix[mask]

    mask = list()
    for st, end in zip(domStartID, domEndID):
        if (np.abs(loopStartID - st) + np.abs(loopEndID - end)).min() < CUTOFF:
            mask.append(False)
        else:
            mask.append(True)
    domainsNoLoops = domains.ix[mask]

    return {'loops': loops, 'domains': domains, 'loopDomains': domainsLoops, 'domainsNoLoops': domainsNoLoops}

def getLoops(mygen, df):
    loops = df
    #loops = pd.read_csv(loops_file, sep="\t")
    loops["chrms"] = [mygen.label2idx[i.strip()[3:]] for i in loops["chr1"].values]
    loops["starts"] = (loops["x1"] + loops["x2"]) / 2
    loops["ends"] = (loops["y1"] + loops["y2"]) / 2
    loops = loops[["chrms", "starts", "ends"]]
    loops = loops.sort_values(["chrms", "starts"])

    return loops


def getDomains(mygen, df):
    domains = df
    #domains = pd.read_csv(domains_file, sep="\t")
    domains["chrms"] = [mygen.label2idx[i[3:]] for i in domains["chr1"].values]
    domains["starts"] = domains["x1"]
    domains["ends"] = domains["x2"]
    domains = domains[["chrms", "starts", "ends"]]
    domains = domains.sort_values(["chrms", "starts"])

    return domains


# load loops, load domains, loop-domains + cache results for retrieval speed
#getLoopsAndDomains = mem.cache(getLoopsAndDomains)
#loopAndDomainDict = {i: getLoopsAndDomains(i) for i in ["hg19", "mm9"]}
#loopAndDomainDict = {i: getLoopsAndDomains(i) for i in ["mm9"]}


# compute average loops
def averageLoops(loopPositions, file_name, pad = 8, doBalance = False):
    c = cooler.Cooler(file_name)
    resolution = c.info['bin-size']
    #pad_bp = pad * resolution
    #pad_bp_m1 = (pad-1) * resolution

    chrmList = list(range(len(c.chromnames)))
    mymaps = list()
    runningCount = 0
    for mychr in chrmList:
        mymap = np.zeros((2 * pad, 2 * pad), np.float64)

        data = c.matrix(balance = doBalance, sparse = True).fetch(c.chromnames[mychr]).tocsr()
        current = loopPositions.ix[loopPositions["chrms"] == mychr]

        if len(current) <= 0:
            #print("len(current) = 0 for chr: {} file: {}".format(mychr, file_name))
            continue

        for (st, end) in zip(current["starts"].values, current["ends"].values):

            if abs(st - end) < 10 * resolution:
                continue

            stBin = int(np.floor(float(st) / resolution))
            endBin = int(np.floor(float(end) / resolution))

            if stBin - pad < 0:
                continue
            if endBin + pad > data.shape[0]:  # len(data):
                continue

            mymap = mymap + data[stBin - pad: stBin + pad, endBin - pad: endBin + pad]
            runningCount += 1

        mymaps.append(mymap)

    return (mymaps, runningCount)

# generate list of shifted loops
def controlLoops(df, numShifted = 100, shiftVal = 500000, excludedRegionSize = 100000):
    dfs = list()
    for i in range(numShifted):
        df = df.copy()
        ran = (np.random.randint(2) * 2 - 1) * (np.random.random(len(df)) * (shiftVal - excludedRegionSize) + excludedRegionSize)
        df["starts"] = df["starts"] + ran
        df["ends"] = df["ends"] + ran
        dfs.append(df)

    df = pd.concat(dfs)
    return df

# generate loop pile-up with controls
def averageLoopsWithControl(loopPositions, filename, numShifted = 5, cg = 1, pad = 8, doBalance = False):
    (mymaps, rc1) = averageLoops(loopPositions, filename, pad = pad, doBalance = doBalance)
    (mymaps2, rc2) = averageLoops(controlLoops(loopPositions, numShifted = numShifted), filename, pad = pad, doBalance = doBalance)

    if not cg == 1:
        mymaps = [coarsegrain(1. * i, cg) / (cg ** 2) for i in mymaps]
        mymaps2 = [coarsegrain(1. * i, cg) / (cg ** 2) for i in mymaps2]

    return (mymaps, mymaps2, rc1, rc2)

# average TAD locations


def averageDomains(file_name, minSize=100000, maxSize=10000000, rescaledDomainSize=30, domainPositions=None):
    # minSize and maxSize limit what "kinds" of domains we're looking for based on size exclusion
    c = cooler.Cooler(file_name)
    resolution = c.info['bin-size']
    chrmList = list(range(len(c.chromnames)))
    maps = list()

    if domainPositions is None:
        assert(domainPositions)

    for mychr in chrmList:
        futureMap = np.zeros(
            (3 * rescaledDomainSize, 3 * rescaledDomainSize), float)

        data = c.matrix(balance=False).fetch(c.chromnames[mychr])

        current = domainPositions.ix[domainPositions["chrms"] == int(mychr)]

        if len(current) <= 0:
            #print("len(current) = 0 for chr: {} file: {}".format(mychr, file_name))
            continue

        for (st, nd) in zip(current["starts"].values, current["ends"].values):
            # exclude by size
            if ((nd - st) < minSize) or ((nd - st) > maxSize):
                continue

            # exclude if too close together (sort of redundant)
            if abs(st - nd) < 10 * resolution:
                continue

            stBin = int(np.rint(float(st) / resolution))
            endBin = int(np.rint(float(nd) / resolution))

            mylen = endBin - stBin + 1
            if stBin - mylen < 0:
                continue

            if endBin + mylen >= len(data):
                continue

            singleMap = 1. * \
                (data[stBin - mylen:endBin + mylen +
                      1, stBin - mylen:endBin + mylen + 1])
            assert(len(singleMap) % 3 == 0)
            try:
                futureMap = futureMap + \
                    zoomArray(singleMap, futureMap.shape, order=1).copy()
            except:
                print('failed at {} {} {} {}'.format(st, nd, mychr, file_name))
                assert(1 == 0)

        maps.append(futureMap)

    return maps
#averageDomains = mem.cache(averageDomains)


# average TAD locations
def averageDomainsError(filename, minSize=100000, maxSize=10000000,
                        rescaledDomainSize=30, domains=None, numShuffles=100):
    # minSize and maxSize limit what "kinds" of domains we're looking for based on size exclusion
    c = cooler.Cooler(filename)
    resolution = c.info['bin-size']
    chrmList = list(range(len(c.chromnames)))
    mymaps = []
    maps = []
    var_maps = []

    if domains is None:
        domains = loopAndDomainDict['mm9']['loops']

    for mychr in chrmList:
        futureMap = np.zeros(
            (3 * rescaledDomainSize, 3 * rescaledDomainSize), float)

        data = c.matrix(balance=False).fetch(c.chromnames[mychr])

        current = domains.ix[domains["chrms"] == mychr]
        # try:
        assert len(current) > 0

        for st, nd in zip(current["starts"].values, current["ends"].values):
            # exclude by size
            if ((nd - st) < minSize) or ((nd - st) > maxSize):
                continue

            stBin, endBin = int(np.rint(float(st) / resolution)
                                ), int(np.rint(float(nd) / resolution))
            mylen = endBin - stBin + 1
            if stBin - mylen < 0:
                continue
            if endBin + mylen >= len(data):
                continue
            singleMap = 1. * \
                (data[stBin - mylen:endBin + mylen +
                      1, stBin - mylen:endBin + mylen + 1])
            assert len(singleMap) % 3 == 0
            futureMap = futureMap + \
                zoomArray(singleMap, futureMap.shape, order=1).copy()
        maps.append(futureMap)

        # do variation
        starts_shuff = current["starts"].values
        ends_shuff = current["ends"].values

        for shuff in range(numShuffles):
            starts_shuff = current["starts"].values
            ends_shuff = current["ends"].values

            rand_vals = np.random.choice(
                np.arange(numShuffles), numShuffles, replace=True)
            starts_shuff = starts_shuff[rand_vals]
            ends_shuff = ends_shuff[rand_vals]

            for st, nd in zip(starts_shuff, ends_shuff):
                # exclude by size
                if ((nd - st) < minSize) or ((nd - st) > maxSize):
                    continue
                stBin, endBin = int(
                    np.rint(float(st) / resolution)), int(np.rint(float(nd) / resolution))
                mylen = endBin - stBin + 1
                if stBin - mylen < 0:
                    continue
                if endBin + mylen >= len(data):
                    continue
                singleMap = 1. * \
                    (data[stBin - mylen:endBin + mylen +
                          1, stBin - mylen:endBin + mylen + 1])
                assert len(singleMap) % 3 == 0
                futureMap = futureMap + \
                    zoomArray(singleMap, futureMap.shape, order=1).copy()
            var_maps.append(futureMap)

        # except:
        #    print("could not process chr{0}".format(mychr))
    return maps, var_maps
#averageDomainsError = mem.cache(averageDomainsError)


# get "GC eigenvector decomposition"
def getGCvector(genome_path, resolution):
    gen = Genome(genome_path, readChrms=["#", "X"])
    gen.setResolution(resolution)
    GC = np.concatenate(gen.GCBin)
    return GC


#getGCvector = mem.cache(getGCvector)

# def chromname2number(chrmname,coolerFile):
#    return [x for x in range(len(coolerFile.chromnames)) if chrmname == coolerFile.chromnames[x]][0]

# def number2chromname(number,coolerFile):
#    return coolerFile.chromnames[number]


def gen2chromname(genLabel, coolerFile):
    # if chromosome contains 'chr', append 'chr'
    if 'chr'+str(genLabel) in coolerFile.chromnames:
        return 'chr'+str(genLabel)
    else:
        return genLabel


def doCoolerSaddle(filename, genome_path, eig = None, gen = 'mm9', desiredResolution = 1e6, matrixSize = 5):
    c = cooler.Cooler(filename)
    file_resolution = getResolution(filename)

    downsampleFactor = int(desiredResolution / file_resolution)
    if downsampleFactor < 1:
        trueResolution = file_resolution

    else:
        trueResolution = file_resolution * downsampleFactor
        if not trueResolution == desiredResolution:
            print("Warning: desiredResolution not possible. Using {} bp resolution instead".format(
                trueResolution))

    # this needs to be changed to work with other genomes
    #genome_path = '/groups/tachibana/sean/data/genomes/mus_musculus/mm9'
    gen = Genome(genome_path, readChrms=["#", "X"])
    gen.setResolution(trueResolution)

    if eig is None:
        eig = np.concatenate(gen.GCBin)

    saddles = list()
    for chrom in range(gen.chrmCount):
        saddle = np.zeros((matrixSize, matrixSize), dtype=float)
        start = gen.chrmStartsBinCont[chrom]
        end = gen.chrmEndsBinCont[chrom]  # - 1

        _chrom = gen2chromname(gen.idx2label[chrom], c)
        cur = c.matrix(balance=False).fetch(_chrom)
        cur = coarsegrain(cur, downsampleFactor, extendEdge=True)
        cur = observedOverExpected(cur)

        mask = np.sum(cur, axis=0) > 0
        cur = cur[mask]
        cur = cur[:, mask]

        GC = eig[start: end]
        GC = GC[mask]

        if len(GC) > matrixSize:
            for i in range(matrixSize):
                for j in range(matrixSize):
                    val = 100//matrixSize
                    (G1, G2) = np.percentile(GC, [val * i, val * i + val])
                    mask1 = (GC > G1) * (GC < G2)
                    (G1, G2) = np.percentile(GC, [val * j, val * j + val])
                    mask2 = (GC > G1) * (GC < G2)
                    saddle[i, j] += cur[np.ix_(mask1, mask2)].mean()
        saddles.append(saddle)

    # print(c.chromnames)
    # sys.exit()
    return saddles
#doCoolerSaddle = mem.cache(doCoolerSaddle)


def doCoolerSaddleError(filename, eig=None, gen='mm9', desiredResolution=1e6, matrixSize=5, numShuffles=100):
    c = cooler.Cooler(filename)

    file_resolution = getResolution(filename)
    downsampleFactor = int(desiredResolution/file_resolution)
    if downsampleFactor < 1:
        trueResolution = file_resolution
    else:
        trueResolution = file_resolution*downsampleFactor
        if trueResolution != desiredResolution:
            print(("Warning: desiredResolution not possible. Using"
                   " {0} bp resolution instead").format(trueResolution))

    if eig is None:
        eig = getGCvector(gen, trueResolution)

    gen = Genome(genome_path, readChrms=["#", "X"])
    gen.setResolution(trueResolution)

    saddles = []
    permuted = []
    saddle = np.zeros((matrixSize, matrixSize), dtype=float)
    for i in range(numShuffles):
        permutted.append(np.zeros((matrixSize, matrixSize), dtype=float))

    for chrom in range(gen.chrmCount):
        st = gen.chrmStartsBinCont[chrom]
        end = gen.chrmEndsBinCont[chrom]
        cur = c.matrix(balance=False).fetch(
            gen2chromname(gen.idx2label[chrom], c))
        cur = coarsegrain(cur, downsampleFactor)
        cur = observedOverExpected(cur)
        mask = np.sum(cur, axis=0) > 0
        cur = cur[mask]
        cur = cur[:, mask]
        GC = eig[st:end]
        GC = GC[mask]
        if len(GC) > matrixSize:
            for i in range(matrixSize):
                for j in range(matrixSize):
                    val = 100//matrixSize
                    G1, G2 = np.percentile(GC, [val * i, val * i + val])
                    mask1 = (GC > G1) * (GC < G2)
                    G1, G2 = np.percentile(GC, [val * j, val * j + val])
                    mask2 = (GC > G1) * (GC < G2)

                    addition = cur[np.ix_(mask1, mask2)]
                    addition = np.reshape(addition, (-1))
                    for k in range(numshuffles):
                        resampled = np.random.choice(
                            addition, len(addition), replace=True)
                        permutted[k][i, j] += resampled.mean()
                    saddle[i, j] += addition.mean()
    return saddle, permutted
#doCoolerSaddleError = mem.cache(doCoolerSaddleError)


def processSaddles(saddles_byChr):
    saddle = functools.reduce(lambda x, y: (x + y), saddles_byChr)
    saddle = saddle / np.mean(saddle)
    saddle = ultracorrect(saddle)
    saddle = np.log2(saddle)
    return saddle


def getMatrixScaling(inMatrix, inMask=list(), measureType='sum', scaleType='log', logFactor=1.3):
    inMatrix = np.array(inMatrix, dtype=np.double)
    N = len(inMatrix)

    if len(inMask) > 0:
        mask2d = inMask
        inMatrix *= mask2d
    else:
        marginals = np.nansum(inMatrix, axis=0)
        mask = marginals > 0
        mask2d = mask[:, None] * mask[None, :]

    if scaleType == 'log':
        bins = logbinsnew(1, N, logFactor)
    else:
        bins = np.arange(0, N)
   
    mids = (0.5 * (bins[:-1] + bins[1:]))
    Pc = list()
    for (st, end) in zip(bins[:-1], bins[1:]):
        curmean = 0
        maskmean = 0
        for i in range(st, end):
            if measureType == 'sum':
                curmean += np.nansum(np.diagonal(inMatrix, i))
                maskmean += np.nansum(np.diagonal(mask2d, i))
            else:
                curmean += np.nanmean(np.diagonal(inMatrix, i))
                maskmean += np.nanmean(np.diagonal(mask2d, i))

        Pc.append(curmean / maskmean)
    mids = np.r_[mids, N]

    #print('len(Pc) = {}'.format(len(Pc)))

    if len(Pc) < 2:
        return (None, None)

    Pc = np.r_[Pc, np.sqrt((Pc[-1] / Pc[-2])) * Pc[-1]]
    return (Pc, mids)


def getCoolerScaling(file_name, reference_cooler_file = None, logFactor = 1.15):
    c = cooler.Cooler(file_name)

    if not reference_cooler_file is None:
        cc = cooler.Cooler(reference_cooler_file)
    else:
        cc = None

    Pc_list = list()
    mids_list = list()
    if not reference_cooler_file is None:
        mapped_chrom = map_chromnames(c.chromnames, cc.chromnames)
    else:
        mapped_chrom = export_chromnames(c.chromnames)

    for (chrom, ref_chrom) in mapped_chrom:
        #print((chrom, ref_chrom))
        # try:
        #inMatrix = c.matrix(balance = False).fetch('{}'.format(chrom))
        inMatrix = c.matrix(balance=False).fetch(chrom)

        # if len(reference_cooler_file) > 0:
        if not reference_cooler_file is None:
            #marginals = np.nansum(cc.matrix(balance = False).fetch('{}'.format(ref_chrom)), axis = 0)
            marginals = np.nansum(
                cc.matrix(balance=False).fetch(ref_chrom), axis=0)
            mask = marginals > 0
            mask2d = mask[:, None] * mask[None, :]
            (Pc, mids) = getMatrixScaling(inMatrix, inMask=mask2d,
                                          measureType='sum', logFactor=logFactor)
        else:
            (Pc, mids) = getMatrixScaling(inMatrix,
                                          measureType='sum', logFactor=logFactor)

        if Pc is None:
            print('Could not process chromosome: {}'.format(chrom))
            continue

        Pc_list.append(Pc)
        mids_list.append(mids)
        # except:
        #print('Could not process chromosome: {}'.format(chrom))

    # get average value genome-wide
    #print('len(Pc_list) = {}'.format(len(Pc_list)))

    biggest_val = np.max([len(x) for x in Pc_list])

    Pc = np.zeros((len(Pc_list), biggest_val)) * np.nan
    for (si, s) in enumerate(Pc_list):
        Pc[si, 0:len(s)] = s

    Pc = np.nanmean(Pc, axis=0)

    mids = mids_list[0]
    for m in mids_list:
        if len(m) > len(mids):
            mids = m

    return (Pc, mids, Pc_list)


def map_chromnames(chromnames, ref_chromnames):
    a = list()

    for ref_chromname in ref_chromnames:
        for chromname in chromnames:
            if chromname[-1 * min(len(chromname), len(ref_chromname)):].lower() == ref_chromname[-1 * min(len(chromname), len(ref_chromname)):].lower():
                a.append((chromname, ref_chromname))
                break

    return a


def export_chromnames(chromnames):
    a = list()

    for chromname in chromnames:
        a.append((chromname, None))

    return a

############### PLOTTING UTILITIES ##########################################


def sumList(matList):
    if len(matList) == 0:
        return 0

    return reduce(lambda x, y: x+y, matList)


def showLoopPileUp(stored_maps, by_keys=None):
    #stored_maps = pickle.load(open(pkl_file, 'rb'))

    if by_keys is None:
        sampleKeys = stored_maps.keys()
    else:
        sampleKeys = by_keys

    for di, dist in enumerate(sampleKeys):
        if di == 0:
            M = sumList(stored_maps[dist][2])
            C = sumList(stored_maps[dist][3])
            m_count = stored_maps[dist][4]
            c_count = stored_maps[dist][5]
        else:
            _M = sumList(stored_maps[dist][2])
            _C = sumList(stored_maps[dist][3])
            _m_count = stored_maps[dist][4]
            _c_count = stored_maps[dist][5]

            M += _M
            C += _C

            if not _m_count is None:
                m_count += _m_count

            if not _c_count is None:
                c_count += _c_count

    return (M, C, m_count, c_count)


def showTADPileUp(stored_tad, downsample_factor=1.0):
    #stored_tad = pickle.load(open(pkl_file, 'rb'))
    for im, m in enumerate(stored_tad):
        if (im == 0):
            M = m
        else:
            M += m
    #diagMaskSize = 3

    if not downsample_factor == 1.0:
        Mshow = np.random.binomial(np.asarray(M, dtype=int), downsample_factor)
    else:
        Mshow = M

    Mshow = np.asarray(Mshow, dtype=float)
    Mshow = (Mshow + Mshow.T) / 2

    ar = np.arange(len(M))
    mask2 = 1 / (1 + np.abs(ar[:, None] - ar[None, :])) ** 0.25

    # do observed over expected
    oe = np.zeros_like(Mshow)
    i, j = np.indices(oe.shape)
    for d in range(len(Mshow)):
        val = np.nanmean(np.diagonal(Mshow, d))
        if val == 0:
            val = 0.01
        oe[i == j-d] = val
        oe[i == j+d] = val
    Mshow_oe = Mshow/oe

    #Mshow_oe = observedOverExpected(Mshow)

    return (Mshow_oe, mask2)


def getNormStrength(M):
    # divide box into 3 equal spaces
    L = len(M)
    box1 = M[0:L//3, 0:L//3]
    box3 = M[L-L//3:L, L-L//3:L]
    return (np.nanmean(box3) + np.nanmean(box1)) / 2


def getLoopStrength(inMatrix, controlMatrix):
    M = inMatrix / controlMatrix
    # divide box into 3 equal spaces
    L = len(M)

    box1 = M[0: L//3, 0: L//3]
    box2 = M[(L//2 - L//6): (L//2 + L//6), (L//2 - L//6): (L//2 + L//6)]
    box3 = M[(L - L//3): L, (L - L//3): L]

    return np.nansum(box2) / (np.nansum(box3) + np.nansum(box1)) * 2


def loopStrap(strengths, numSamples=10000, conf_level=0.95):
    N = len(strengths)

    # randomly sample strengths numSamples times, and calculate loop strength
    loop_base = np.zeros(numSamples)
    for i in range(numSamples):
        randN = np.random.choice(N, N, replace=True)
        M = reduce((lambda x, y: x + y), [strengths[x][0] for x in randN])
        C = reduce((lambda x, y: x + y), [strengths[x][1] for x in randN])
        loop_base[i] = getLoopStrength(M, C)  # ,pad=pad)

    s_loop = np.sort(loop_base)
    lo = int(numSamples*(1-conf_level))
    hi = int(numSamples*(conf_level))
    med = int(numSamples*0.5)
    ci_lower = s_loop[med] - s_loop[lo]
    ci_higher = s_loop[hi] - s_loop[med]

    # get true loop strength
    randN = np.random.choice(N, N, replace=False)
    M_true = reduce((lambda x, y: x + y), [strengths[x][0] for x in randN])
    C_true = reduce((lambda x, y: x + y), [strengths[x][1] for x in randN])

    return getLoopStrength(M_true, C_true), np.nanstd(loop_base), ci_lower, ci_higher


def _getCompartmentStrength(M):
    L = len(M)-1
    return np.log(M[0, 0]*M[L, L] / (M[0, L]*M[L, 0]))


def _compartmentStrap(saddles, numSamples=10000, conf_level=0.95):
    N = len(saddles)
    # randomly sample strengths numSamples times, and calculate loop strength
    cs = np.zeros(numSamples)
    for i in range(numSamples):
        randN = np.random.choice(N, N, replace=True)
        M = reduce((lambda x, y: x + y), [saddles[x] for x in randN])
        M = M / np.mean(M)
        M = ultracorrect(M)
        cs[i] = getCompartmentStrength(M)  # ,pad=pad)

    s_cs = np.sort(cs)
    lo = int(numSamples*(1-conf_level))
    hi = int(numSamples*(conf_level))
    med = int(numSamples*0.5)
    ci_lower = s_cs[med] - s_cs[lo]
    ci_higher = s_cs[hi] - s_cs[med]

    # get true loop strength
    randN = np.random.choice(N, N, replace=False)
    M_true = reduce((lambda x, y: x + y), [saddles[x] for x in randN])
    M_true = M / np.mean(M_true)
    M_true = ultracorrect(M_true)

    return getCompartmentStrength(M_true), np.nanstd(cs), ci_lower, ci_higher


def getTADStrength(M):
    L = len(M)
    box1 = 0.5*np.sum(M[0:L//3, L//3:2*L//3])+0.5 * \
        np.sum(M[L//3:2*L//3, 2*L//3:L])
    box2 = np.sum(M[L//3:2*L//3, L//3:2*L//3])
    return box2 / box1


############### UTILITY TOOLS FOR MERGING COOLER FILES ######################


class _CoolerMerger():
    def __init__(self, listOfCoolers, **kwargs):
        self.c_list = listOfCoolers
        self.matrix_list = list()

        for c in listOfCoolers:
            temp = c.pixels()[:]
            self.matrix_list.append(temp)

    def aggregate(self):
        table = pd.concat(self.matrix_list, ignore_index=True)
        return (table.groupby(['bin1_id', 'bin2_id']))['count'].sum().reset_index()

    def size(self):
        return np.sum(c.info['nnz'] for c in self.c_list)

    def __iter__(self):
        df = self.aggregate()
        yield {k: v.values for (k, v) in six.iteritems(df)}


def _collectMergedCoolerMetaData(listOfCoolers):
    new_metadata = []
    for ic, c in enumerate(listOfCoolers):
        if ic == 0:
            new_metadata = {k: int(c.info['metadata'][k])
                            for k in c.info['metadata'].keys()}
        else:
            new_metadata = {k: (int(
                c.info['metadata'][k]) + int(new_metadata.get(k, 0))) for k in c.info['metadata'].keys()}
    if 'files_merged' not in new_metadata.keys():
        new_metadata['files_merged'] = len(listOfCoolers)
    return new_metadata


def mergeCoolerFiles(data_files_list, merged_filename):
    listOfCoolers = list()
    for data_file in data_files_list:
        try:
            listOfCoolers.append(cooler.Cooler(data_file))
        except:
            raise ValueError

    iterator = _CoolerMerger(listOfCoolers)
    c = listOfCoolers[0]

    print('merging: {}'.format(os.path.split(merged_filename)[1]))
    # cooler.io.create(merged_filename, c.bins()[:], pixels = iterator, metadata=_collectMergedCoolerMetaData(
    #    listOfCoolers), assembly=c.info['genome-assembly'])

    cooler.create_cooler(merged_filename, c.bins()[:], pixels=iterator, metadata=_collectMergedCoolerMetaData(
        listOfCoolers), assembly=c.info['genome-assembly'])


def natural_sort(l):
    def convert(text): return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): return [convert(c)
                                   for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)
# http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
