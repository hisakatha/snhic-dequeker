#! /usr/bin/env python
# (c) 2019. All Rights Reserved.
# Code written by: Hugo Brandao (hbrandao@g.harvard.edu)
# Code written by: Maksim Imakaev (imakaev@mit.edu)
# Code written by: Sean Powell (sean.powell@imba.oeaw.ac.at)

# NOTE: this version of the code only downsamples analyses for:
# the loop, TAD, saddle strengths
# it currently does not downsample the analysis for scaling plots


from lib import defaults
__version__ = defaults.VERSION

import os
import sys
import click
import pickle
import datetime

import cooler
import matplotlib
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import seaborn as sns

from lib import coolerAnalysis as css
from lib import snhicAnalysis as ass
from lib.common import read_samples_file, get_chunks, group_by_tags, get_cooler_files
from mirnylib.numutils import observedOverExpected, ultracorrect

from mirnylib import numutils


from mpl_toolkits.axes_grid1 import make_axes_locatable
from pandas import ExcelWriter
from scipy.ndimage.filters import gaussian_filter1d
from cooltools import saddle

matplotlib.use('Agg')
DOWNSAMP_READ_COUNTS = -1 # defaults to minimum of data sets
SMOOTH_SIGMA = defaults.SMOOTH_SIGMA
RESOLUTION = defaults.RESOLUTION
SAVE_SOURCE_DATA = defaults.SAVE_SOURCE_DATA
MIN_EXP_LOOP_SIZE = defaults.MIN_EXP_LOOP_SIZE
MAX_EXP_LOOP_SIZE = defaults.MAX_EXP_LOOP_SIZE
METAKEYS = defaults.METAKEYS
INCLUDE_COLUMN = defaults.INCLUDE_COLUMN

COOL_EXT = defaults.COOLER_EXT
OUTPUT_EXT = '.{resolution}{output_ext}'

TAGS = defaults.TAG
PLOTS_PER_PAGE = defaults.PLOTS_PER_PAGE

SADDLE_PKL_FILE = defaults.SADDLE_PKL_FILE
PILEUP_LOOPS_PKL_FILE = defaults.PILEUP_LOOPS_PKL_FILE
SCALINGCOOLER_PKL_FILE = defaults.SCALINGCOOLER_PKL_FILE
REGIONAVG_DOMAINS_PKL_FILE = defaults.REGIONAVG_DOMAINS_PKL_FILE

PDF_FILE = '{}{}'
STATS_FILE = '{}{}'

@click.version_option(version = __version__)
@click.command(context_settings = defaults.CONTEXT_SETTINGS)
@click.option('--samples-file',
               type = str,
               required = True)
@click.option('--data-processed-dir',
               type = str,
               required = True)
@click.option('--output-dir',
               type = str,
               required = True)
@click.option('--tags',
               type = str,
               default = TAGS,
               show_default = True,
               required = True)
@click.option('--stats-dir',
               required = True,
               type = str,
               help = 'The directory for the statistics.')
@click.option('--cooler-dir',
               required = True,
               type = str,
               help = 'The directory with the merged cooler files.')
@click.option('--cooler-ext',
               default = COOL_EXT,
               show_default = True,
               type = str,
               help = 'The cooler file extension.')
@click.option('--downsamp-readcount',
               default = DOWNSAMP_READ_COUNTS,
               show_default = True,
               type = int,
               help = 'The number of reads to use for downsampling.')
@click.option('--figure-name',
               type = str,
               default = 'Figure_1',
               show_default = True,
               help = 'The name of the figure.')
@click.option('--resolution',
               type = int,
               default = RESOLUTION,
               show_default = True,
               help = 'The resolution of the cooler files that will be analzed.')
@click.option('--min-exp-loop-size',
               type = int,
               default = MIN_EXP_LOOP_SIZE,
               show_default = True,
               help = 'The minimum expected loop size.')
@click.option('--max-exp-loop-size',
               type = int,
               default = MAX_EXP_LOOP_SIZE,
               show_default = True,
               help = 'The maximum expected loop size.')
@click.option('--plots-per-page',
               type = int,
               default = PLOTS_PER_PAGE,
               show_default = True,
               help = 'The amount of samples per page page.')
@click.option('--include-column',
               default = INCLUDE_COLUMN,
               show_default = True,
               type = str,
               help = 'The column in the samples file that determines which samples to include in the analyses.')
@click.option('--balance/--no-balance',
               is_flag = True,
               help = 'Set if the contact matrices were balanced. Default: not set')
@click.option('--balance-saddle',
               is_flag = True,
               help = 'Set if the contact matrices were not balanced, but you want the final saddle plot to be balanced. This has no effect if --saddle-plots is not set and/or --balance is set. Default: not set')
@click.option('--saddle-plots/--no-saddle-plots',
               is_flag = True,
               help = 'Create the saddle plots. Default: --no-saddle-plots')
@click.option('--by-chrom/--by-genome',
               is_flag = True,
               help = 'Analyse by chromosome or the whole genome. Default: --by-genome')
def main(samples_file, data_processed_dir, output_dir, tags, stats_dir, cooler_dir, cooler_ext, figure_name, resolution, min_exp_loop_size, max_exp_loop_size, plots_per_page, include_column, balance, balance_saddle, saddle_plots, by_chrom, downsamp_readcount):
    '''Create the summary PDFs.
    '''

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    assert(os.path.exists(output_dir))


    df = read_samples_file(samples_file)
    df = df[df[include_column] == True]

    if len(df) == 0:
        sys.stderr.write('error: --samples-file argument of the incorrect file type. Only support for .csv, .xls, .xlsx.')
        sys.exit()

    tag_names = tags.split(',')
    grouped_dict = group_by_tags(df, tag_names)

    merged_cooler_files = get_cooler_files(grouped_dict, resolution, cooler_dir, cooler_ext)

    ################# AUTOMATIC DOWNSAMPLING FACTOR CALCULATION ################
    # get downsampling probabilities for all cooler files
    cis_count_list = []
    for merged_cooler_file in merged_cooler_files:
        c = cooler.Cooler(merged_cooler_file)
        meta_dict = c.info['metadata']
        cis_contacts = meta_dict.get('cis', 0)
        cis_count_list.append(cis_contacts)

    if downsamp_readcount < 0:
        downsampled_cis_readcount = np.min(cis_count_list)
    else:
        downsampled_cis_readcount = desired_cis_readcount

    p_binom_downsampling = [downsampled_cis_readcount/x for x in cis_count_list]
    ################# END DOWNSAMPLING FACTOR CALCULATION #################

    if SAVE_SOURCE_DATA:
        xls_path = os.path.join(output_dir, 'source_data_{}.xls'.format(figure_name))
        writer = ExcelWriter(xls_path)
    else:
        writer = None

    stats_file = os.path.join(stats_dir, STATS_FILE.format(figure_name, defaults.FINAL_STATS_EXT))
    pdf_name = os.path.join(output_dir, PDF_FILE.format(figure_name, defaults.PDF_EXT))

    metadatas = list()
    dimX, dimY = (0, 6)
    with PdfPages(pdf_name) as pdf:
        for merged_cooler_page_files in get_chunks(merged_cooler_files, plots_per_page):
            (_, pdf_file_name) = os.path.split(pdf_name)

            num_samples = len(merged_cooler_page_files)
            dimX, dimY = (max(dimX, num_samples), dimY)
            fig = plt.figure(figsize = (3 * dimX, 3 * dimY))
            plt.suptitle(pdf_file_name, fontsize = 8)
            panels = gs.GridSpec(dimY, dimX)

            for (si, merged_cooler_file) in enumerate(merged_cooler_page_files):
                (_, merged_cooler_file_name) = os.path.split(merged_cooler_file)
                sheet_name = merged_cooler_file_name.replace('-', '_')
                sheet_name = merged_cooler_file_name.split('.')[0]

                (_, merged_cooler_file_name) = os.path.split(merged_cooler_file)

                plot_title = merged_cooler_file_name.split('.')[0]
                print('no1')
                row_num = 0
                pileup_loops_pkl_file_name = PILEUP_LOOPS_PKL_FILE.format(merged_cooler_file_name)
                pileup_loops_pkl_file = os.path.join(data_processed_dir, pileup_loops_pkl_file_name)
                stored_loops = pickle.load(open(pileup_loops_pkl_file, 'rb'))
                loop_strength = doLoops_downsampled(plt, panels, plot_title, stored_loops, num_samples, resolution, writer, sheet_name, row_num, si, p_binom_downsampling)

                print('no2')
                row_num += 1
                regionalavg_domains_pkl_file_name = REGIONAVG_DOMAINS_PKL_FILE.format(merged_cooler_file_name)
                regionalavg_domains_pkl_file = os.path.join(data_processed_dir, regionalavg_domains_pkl_file_name)
                stored_tads = pickle.load(open(regionalavg_domains_pkl_file, 'rb'))
                tad_strength = doTADS(plt, panels, stored_tads, num_samples, writer, sheet_name, row_num, si,p_binom_downsampling)

                print('no3')
                saddle_strength = None
                print(numutils.__file__)
                if saddle_plots:
                    row_num += 1
                    saddle_pkl_file_name = SADDLE_PKL_FILE.format(merged_cooler_file_name)
                    saddle_pkl_file = os.path.join(data_processed_dir, saddle_pkl_file_name)
                    stored_saddle = pickle.load(open(saddle_pkl_file, 'rb'))
                    saddle_strength = doSaddles(plt, panels, stored_saddle, num_samples, writer, sheet_name, row_num, si,p_binom_downsampling)

                #print('no4')
                row_num += 1
                scaling_pkl_file_name = SCALINGCOOLER_PKL_FILE.format(merged_cooler_file_name)
                scaling_pkl_file = os.path.join(data_processed_dir, scaling_pkl_file_name)
                stored_scaling = pickle.load(open(scaling_pkl_file, 'rb'))
                vals = doScaling(plt, panels, stored_scaling, SMOOTH_SIGMA, resolution, row_num, si)

                #print('no5')
                row_num += 1
                (loop_size_from, loop_size_to) = doSlopes(plt, panels, vals, SMOOTH_SIGMA, resolution, writer, sheet_name, row_num, si, min_exp_loop_size, max_exp_loop_size)
                #print('no6')
                row_num += 1
                metadata = doStatsAndMetadata(plt, panels, merged_cooler_file, merged_cooler_file, loop_size_from, loop_size_to, resolution, writer, sheet_name, row_num, si,p_binom_downsampling)
                print('no7')

                metadata['id'] = plot_title
                metadata['loop_strength'] = loop_strength
                metadata['tad_strength'] = tad_strength
                if saddle_plots:
                    metadata['saddle_strength'] = saddle_strength

                metadatas.append(metadata)


            plt.draw()
            pdf.savefig()
            plt.close()

        first_columns = ['id', 'files_merged', 'loop_strength', 'tad_strength']
        if saddle_plots:
            first_columns.append('saddle_strength')

        first_columns.append('avg_loop_size')

        export_metadata(metadatas, stats_file, first_columns)

        d = pdf.infodict()
        d['Author'] = 'Mirny and Tachibana Lab pairtools and snHi-C pipeline'
        d['Keywords'] = 'Hi-C, snHi-C'
        d['CreationDate'] = datetime.datetime.today()

    if SAVE_SOURCE_DATA:
        writer.save()

def export_metadata(metadatas, full_file_name, first_columns):
    columns = list()
    columns.extend(first_columns)
    for column_name in metadatas[0].keys():
        if column_name in first_columns:
            continue

        columns.append(column_name)

    arr = list()
    for metadata in metadatas:
        a = list()
        for column in columns:
            a.append(metadata.get(column, 'NA'))

        arr.append('\t'.join([str(i) for i in a]))


    f_out = open(full_file_name, 'w')
    f_out.write('\t'.join(columns))
    f_out.write('\n')
    f_out.write('\n'.join(arr))
    f_out.close()

def plot_saddles(fig, plt, panel, clr, balance_saddle, long_name, resolution, data_folder, panel_num, si, histbins = 6, quantile_binning = True, balance = False):
     #######################################################################
    # the "saddle" plot
    FONTSIZE = 8

    #print('saddle plot {}'.format(long_name))
    exp_file_name = '{name}.{res}.expected.cis.tsv'.format(name = long_name, res = resolution)
    exp_file = os.path.join(data_folder, exp_file_name)
    eig_file_name = '{name}.{res}.eigs.cis.vecs.txt'.format(name = long_name, res = resolution)
    eig_file = os.path.join(data_folder, eig_file_name)

    exp = pd.read_table(exp_file)
    eig = pd.read_table(eig_file)

    # Determine how to bin the range of the E1 signal
    if quantile_binning:
        q_binedges = np.linspace(0, 1, histbins)
        binedges = saddle.quantile(eig['E1'], q_binedges)
    else:
        (qlo, qhi) = saddle.quantile(eig['E1'], [0.02, 0.98])  # trim outliers
        binedges = np.linspace(qlo, qhi, histbins)

    # Digitize the signal into integers
    (digitized, hist) = saddle.digitize_track(binedges, track = (eig, 'E1'))

    # Construct a function that fetches and calculates observed/expected
    #print('observed / expected')
    if balance:
        getmatrix = saddle.make_cis_obsexp_fetcher(clr, (exp, 'balance.avg'))
    else:
        getmatrix = ass.make_cis_obsexp_fetcher(clr, (exp, 'count.avg'))

    # Build the saddle histogram
    #print('saddle histogram')
    pkl_file_name = '{name}.{res}.saddledata.pkl'.format(name = long_name, res = resolution)
    pkl_file = os.path.join(data_folder, pkl_file_name)
    if not os.path.exists(pkl_file):
        (sums, counts) = saddle.make_saddle(getmatrix, binedges, (digitized, 'E1.d'), contact_type = 'cis')
        saddle_data = sums / counts
        f = open(pkl_file, 'wb')
        #print('writing: {}'.format(f.name))
        pickle.dump((sums, counts), f)

    else:
        f = open(pkl_file, 'rb')
        #print('loading: {}'.format(f.name))
        (sums, counts) = pickle.load(f)
        saddle_data = sums / counts

    #print('saddle_data:')
    #print(saddle_data)
    if balance_saddle:
        saddle_data = ass.balance_saddle_data(saddle_data)

    AA = saddle_data[1][1]
    AB = saddle_data[1][histbins - 1]
    BA = saddle_data[histbins - 1][1]
    BB = saddle_data[histbins - 1][histbins - 1]

    #saddle_data = np.log2(saddle_data)
    #print(saddle_data)

    cs = np.log10((AA * BB ) / (AB * BA))
    cs = (AA * BB ) / (AB * BA)

    pal = sns.color_palette('Greys')
    # Make the saddle plot
    #print('generating plot')
    if quantile_binning:
        if si == None:
            g = saddle.saddleplot(q_binedges, hist, saddle_data, heatmap_kws = {'vmin': -1.0, 'vmax': 1.0}, color = pal[0], fig = fig, subplot_spec = panel[panel_num])
        else:
            g = saddle.saddleplot(q_binedges, hist, saddle_data, heatmap_kws = {'vmin': -1.0, 'vmax': 1.0}, color = pal[0], fig = fig, subplot_spec = panel[panel_num, si])
    else:
        if si == None:
            g = saddle.saddleplot(binedges, hist, saddle_data, heatmap_kws = {'vmin': -1.0, 'vmax': 1.0}, color = pal[0], fig = fig, subplot_spec = panel[panel_num])
        else:
            g = saddle.saddleplot(binedges, hist, saddle_data, heatmap_kws = {'vmin': -1.0, 'vmax': 1.0}, color = pal[0], fig = fig, subplot_spec = panel[panel_num, si])


    offsetX = 0.05
    offsetY = 0.025
    offsetMidX = 0.03
    #plt.title(long_name, loc = 'left')
    g['ax_heatmap'].text(0.1 - offsetX , 0.1 + offsetY, "AA", fontsize = FONTSIZE)
    g['ax_heatmap'].text(0.9 - offsetX , 0.1 + offsetY, "AB", fontsize = FONTSIZE)
    g['ax_heatmap'].text(0.1 - offsetX , 0.9 + offsetY, "BA", fontsize = FONTSIZE)
    g['ax_heatmap'].text(0.9 - offsetX , 0.9 + offsetY, "BB", fontsize = FONTSIZE)
    g['ax_heatmap'].text(0.5 - offsetX - offsetMidX , 0.5 + offsetY, '{:.2f}'.format(cs), fontsize = FONTSIZE)
    g['ax_margin_y'].set_visible(False)
    g['ax_margin_x'].set_visible(False)

    return (sums, counts)

def doLoops_downsampled(plt, panels, plot_title, stored_maps, num_samples, resolution, writer, sheet_name, row_num, si, p_binom_downsampling):
    # DO LOOPS
    plt.subplot(panels[row_num, si]) #  (pi+1)*(si+1)
    plt.title(plot_title, fontsize = 10)

    colorMax = 1
    (M, C, m_count, c_count) = css.showLoopPileUp(stored_maps)

    # code added for downsampling
    p_downsample = p_binom_downsampling[si]
    M = np.random.binomial(np.asarray(M,dtype=int),p_downsample)
    C = np.random.binomial(np.asarray(C,dtype=int),p_downsample)

    if m_count == 0 or c_count == 0:
        strength = 0
    else:
        strength = css.getLoopStrength(M / m_count, C / c_count)

    if m_count == 0:
        Mshow = 0
    else:
        Mshow = M / (C)  * c_count / m_count

    Mshow = Mshow / css.getNormStrength(Mshow)

    im = plt.imshow(np.log2(Mshow), cmap = 'coolwarm', vmax = colorMax, vmin = -colorMax)
    dist_lo = int(resolution / 1e3 * len(M) / 2)
    dist_hi = int(resolution / 1e3 * ((len(M) / 2) - 1))

    ax = plt.gca()
    ax.annotate(str(round(strength, 3)), xy =(2 * len(Mshow) / 3, 2 * len(Mshow) / 3), xytext = (2 * len(Mshow) / 3, 2 * len(Mshow) / 3))

    if (si == 0) :
        plt.ylabel('Loops')
        plt.yticks([0, len(Mshow)//2, len(Mshow)-1], ["-{0}kb".format(str(dist_lo)), "0 kb", "+{0}kb".format(str(dist_hi))], fontsize = 8)
        plt.xticks([0, len(Mshow)//2, len(Mshow)-1], ["-{0}kb".format(str(dist_lo)), " ", "+{0}kb".format(str(dist_hi))], fontsize = 8)
    else:
        plt.yticks([0, len(Mshow)//2, len(Mshow) - 1], ["", "", ""], fontsize = 8)
        plt.xticks([0, len(Mshow)//2, len(Mshow) - 1], ["", "", ""], fontsize = 8)

    if si == num_samples - 1:
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size = "5%", pad = 0.1)
        cbar = plt.colorbar(im, cax = cax, ticks = [-colorMax, colorMax])
        cbar.set_label('log2 (enrichment)', rotation = 90)

    if SAVE_SOURCE_DATA:
        Mshow_df = pd.DataFrame(Mshow)
        Mshow_df.to_excel(writer,'Loops_' + sheet_name[0:25])

    return strength

def doTADS(plt, panels, stored_maps, num_samples, writer, sheet_name, row_num, si,p_binom_downsampling):
    # DO TADS
    vmin = 0.3
    vmax = 0.9

    plt.subplot(panels[row_num, si]) #  (pi+1)*(si+1)
    colorMax = 1

    (Mshow_oe, scaling) = css.showTADPileUp(stored_maps,p_binom_downsampling[si])

    tadStrength = css.getTADStrength( observedOverExpected(Mshow_oe * scaling) )

    im = plt.imshow((Mshow_oe * scaling), interpolation = 'nearest', cmap = 'hot_r', vmin = vmin, vmax = vmax)

    if SAVE_SOURCE_DATA:
        Mshow_df = pd.DataFrame(Mshow_oe * scaling)
        Mshow_df.to_excel(writer, 'TADs_'+ sheet_name[0: 25])

    if  si == 0:
        plt.yticks([len(Mshow_oe) / 3, len(Mshow_oe) / 2, 2 * len(Mshow_oe) / 3], ["", "TAD", ""], rotation = 'vertical', va = 'center')
        plt.xticks([len(Mshow_oe) / 3, len(Mshow_oe) / 2, 2 * len(Mshow_oe) / 3], ["", "TAD", ""], rotation = 'horizontal',va = 'center')
    else:
        plt.yticks([len(Mshow_oe) / 3, 2 * len(Mshow_oe) / 3], ["", ""])
        plt.xticks([len(Mshow_oe) / 3, 2 * len(Mshow_oe) / 3], ["", ""])

    plt.annotate(str(round(tadStrength, 3) ),xy = (2 * len(Mshow_oe) / 3, len(Mshow_oe) / 2), xytext = (2 * len(Mshow_oe) / 3, len(Mshow_oe) / 2))

    if si == 0:
        plt.ylabel('Domains')

    if si == (num_samples - 1):
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size = "5%", pad = 0.05)
        cbar = plt.colorbar(im, cax = cax, ticks = [vmin, vmax])
        cbar.set_label('contact enrichment', rotation = 90)

    return tadStrength

def doSaddles(plt, panels, stored_saddles, num_samples, writer, sheet_name, row_num, si,p_binom_downsampling):
    # DO SADDLES
    ####################
    # need to integrate this row_num
    plt.subplot(panels[row_num, si]) #  (pi+1)*(si+1)

    colorMax = 0.8
    M = list()

    #stored_saddles = pickle.load(open(saddle_pkl_file, 'rb'))

    M = css.sumList(stored_saddles)

    # apply downsampling
    M = np.random.binomial(np.asarray(M,dtype=int),p_binom_downsampling[si])
    M = (M+M.T)/2

    M = M / np.mean(M)
    M = ultracorrect(M)

    cs = np.log(M[0, 0] * M[4, 4] / (M[0, 4] * M[4, 0]))

    M = np.log2(M)
    im = plt.imshow(M, interpolation = 'nearest', cmap = 'coolwarm', vmin = -colorMax, vmax = colorMax)

    if SAVE_SOURCE_DATA:
        Mshow_df = pd.DataFrame(M)
        Mshow_df.to_excel(writer, 'Comp_' + sheet_name[-24:])

    if si == 0:
        plt.ylabel("Compartments")
        plt.annotate("AA", xy = (-0.25, 0.25))
        plt.annotate("AB", xy = (4 - 0.25, 0.25))
        plt.annotate("BA", xy = (-0.3, 4 + 0.25))
        plt.annotate("BB", xy = (4 - 0.3, 4 + 0.25))
    plt.annotate( str(round(cs, 3)), xy = ((len(M) / 2) - 1, len(M) / 2), xytext = ((len(M) / 2) - 1, (len(M) / 2) - 0.25))

    if  si == 0:
        plt.yticks([0, len(M) - 1], ["active", "inactive"],rotation = 'vertical', va = 'center')
        plt.xticks([0, len(M) - 1], ["active", "inactive"])
    else:
        plt.yticks([0, len(M) - 1], ["", ""], rotation = 'vertical', va = 'center')
        plt.xticks([0, len(M) - 1], ["", ""])

    if si == (num_samples - 1):
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size = "5%", pad = 0.05)
        cbar = plt.colorbar(im, cax = cax, ticks=[-colorMax, colorMax])
        cbar.set_label('log2 (enrichment)', rotation = 90)

    return cs

def doScaling(plt, panels, vals, smooth_sigma, resolution, row_num, si):
    # DO SCALING

    plt.subplot(panels[row_num, si]) #  (pi+1)*(si+1)

    #vals = pickle.load(open(scaling_pkl_file, 'rb'))

    plt.loglog(vals[1] * resolution, vals[0] / vals[0][1])
    plt.semilogx(gaussian_filter1d(vals[1] * resolution, smooth_sigma), gaussian_filter1d( vals[0] / vals[0][1], smooth_sigma ))

    plt.gca().set_aspect('equal', adjustable = 'box')
    plt.grid()
    plt.ylim([1e-5, 1.5])
    plt.xlim([1e4, 1e8])
    if si == 0:
        plt.ylabel('Contact probability, $P_c(s)$')
        plt.xlabel('Genomic separation (bp)')
    else:
        plt.gca().yaxis.set_ticklabels([])

    return vals


def doSlopes(plt, panels, vals, smooth_sigma, resolution, writer, sheet_name, row_num, si, min_exp_loop_size, max_exp_loop_size):
    # DO SLOPES

    # vals[0] = Pc or Pc_list
    # vals[1] = mids

    plt.subplot(panels[row_num, si]) #  (pi+1)*(si+1)
    dlogx = np.diff(np.log(vals[1]))
    slope = np.diff(np.log(vals[0])) / dlogx
    plt.semilogx(vals[1][1:] * resolution, slope)
    plt.semilogx(gaussian_filter1d(vals[1][1:] * resolution, smooth_sigma), gaussian_filter1d(slope, smooth_sigma))

    #loop_size_filt_idx = np.argmax(gaussian_filter1d(slope[(vals[1][1:] * resolution < max_exp_loop_size)], smooth_sigma))

    mask = (vals[1][1:] * resolution > min_exp_loop_size) & (vals[1][1:] * resolution < max_exp_loop_size)
    mask = np.logical_not(mask)

    _slope = np.copy(slope)
    _slope = gaussian_filter1d(_slope, smooth_sigma)

    #print(_slope)
    _slope[mask] = -1000

    loop_size_filt_idx = np.argmax(_slope)
    loop_size_from = vals[1][loop_size_filt_idx - 1]
    loop_size_to = vals[1][loop_size_filt_idx]

    #print(loop_size_filt_idx)

    #print('{:.0f}-{:.0f} kbp'.format(loop_size_from * resolution / 1e3, loop_size_to * resolution / 1e3))

    plt.gca().set_aspect('equal', adjustable = 'box')
    plt.grid()
    plt.ylim([-3,0])
    plt.xlim([1e4,1e8])

    if si == 0:
        plt.ylabel('Slope, log($P_c(s)$)/log($\delta$ s)')
        plt.xlabel('Genomic separation (bp)')

    if SAVE_SOURCE_DATA:
        scaling_df = pd.DataFrame()
        scaling_df['Genomic distance (bp)'] = vals[1] * resolution
        scaling_df['Contact probability '] = vals[0] / vals[0][1]
        scaling_df['log(slope)'] = np.r_[0, slope]
        scaling_df.to_excel(writer, 'Pc(s)_' + sheet_name[0:25])

    return (loop_size_from, loop_size_to)

def doStatsAndMetadata(plt, panels, sample_keys, cooler_file, loop_size_from, loop_size_to, resolution, writer, sheet_name, row_num, si,p_binom_downsampling):
    # DO SHOW SUMMARY STATISTICS + METADATA

    #plt.subplot(panels[row_num + 1, si]) #  (pi+1)*(si+1)
    plt.subplot(panels[row_num, si]) #  (pi+1)*(si+1)


    c = cooler.Cooler(cooler_file)

    # add loop size statistic to file
    meta_dict = c.info['metadata']
    meta_dict['avg_loop_size'] = '{:.0f}-{:.0f} kbp'.format(loop_size_from * resolution / 1e3, loop_size_to * resolution / 1e3)
    disp_text = list()

    cis_contacts = int(meta_dict.get('cis', 0))
    total_contacts = int(meta_dict.get('total', 0))


    meta_dict['cis_post_downsampling'] = int(cis_contacts*p_binom_downsampling[si])
    meta_dict['downsampling_factor'] = p_binom_downsampling[si]

    #print(isinstance(meta_dict['total'], int))
    for i in meta_dict:

        if any([x == i for x in METAKEYS]):
            try:
                disp_text.append((i, '{:,}'.format(int(meta_dict[i]))))#text(i, y, s, fontsize=12)
            except:
                disp_text.append((i, '{}'.format(meta_dict[i])))

            if 'cis_' in i and cis_contacts > 0:
                disp_text.append((i, '{:.2f}%'.format(100 * meta_dict[i]*p_binom_downsampling[si] / cis_contacts)))

        if i == 'downsampling_factor' and total_contacts > 0:
            disp_text.append((i, '{:.2f}'.format(meta_dict[i])))

        if i == 'cis_post_downsampling' and total_contacts > 0:
            disp_text.append((i, '{:.2f}'.format(int(meta_dict[i]))))

        if i == 'total_nodups' and total_contacts > 0:
            disp_text.append((i, '{:.2f}%'.format(100 * int(meta_dict[i]*p_binom_downsampling[si]) / total_contacts)))

        if i == 'cis' and total_contacts > 0:
            disp_text.append((i, '{:.2f}%'.format(100 * int(meta_dict[i]*p_binom_downsampling[si]) / total_contacts)))

        if i == 'trans' and total_contacts > 0:
            disp_text.append((i, '{:.2f}%'.format(100 * int(meta_dict[i]*p_binom_downsampling[si]) / total_contacts)))



    plt.text(0,0,"\n".join(["{} = {}".format(x[0], x[1]) for x in disp_text]), fontsize = 8)
    plt.axis('off')

    if SAVE_SOURCE_DATA:
        stats_df = pd.DataFrame.from_dict(meta_dict, orient = 'index')
        stats_df.to_excel(writer, 'Stats_' + sheet_name[0:25])

    return meta_dict

if __name__ == "__main__":
    main()
