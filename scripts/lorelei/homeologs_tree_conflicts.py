#!/usr/bin/env python

"""
    Barplots and hypergeomtric tests.

    Example::

        $ python -m scripts.lorelei.homeologs_tree_conflicts TODO
"""
import sys
import argparse

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as ss

from statsmodels.sandbox.stats.multicomp import multipletests


def load_counts(input_file):

    """
    Loads input data for hypergeom test: each line is a count category pair, space-separated.

    Args:
        input_file (str) : path to the input file

    Returns:
        dict: for each category (chromosomes), as keys, the count (value)
    """

    data = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            if line[0] != "#":
                counts, chrom = line.strip().split()
                try:
                    chrom = int(chrom)
                except ValueError:
                    pass
                data[chrom] = int(counts)
    return data


def hypergeom_enrich_depl(data_obs, data_tot, alpha=0.05, multitest_adjust='fdr_bh'):

    """
    Hypergeometric tests for enrichment and depletion with multiple testing correction.

    Args:
        data_obs (dict): for each category, observed counts
        data_tot (dict): for each category, total number of objects
        alpha (float): significance level
        multitest_adjust (str): method to adjust pvalues for multiple testing

    Returns:
        tuple of lists: categories, corresponding proportion of observed counts, enrichment or
        depletion, whether null hyp. is rejected, and adjusted p-values.
    """

    values = []
    direction = []
    pvalues = []
    tot = sum(data_tot.values())
    tot_obs = sum(data_obs.values())
    for chrom in sorted(list(data_obs.keys())):

        exp = (tot_obs/tot)*data_tot[chrom]
        obs = data_obs[chrom]
        if obs > exp:
            pval = ss.hypergeom.sf(obs-1, tot, tot_obs, data_tot[chrom]) #enrichment
            direction.append('enrichment')
        else:
            pval = ss.hypergeom.cdf(obs+1, tot, tot_obs, data_tot[chrom]) #depletion
            direction.append('depletion')

        pvalues.append(pval)

        values.append(data_obs[chrom]/data_tot[chrom]*100)

    mtest = multipletests(pvalues, alpha=alpha, method=multitest_adjust)

    return (sorted(list(data_obs.keys())), values, direction, list(mtest[0]), list(mtest[1]))


def plot_sign(ax, to_highlight=None):
    """
    Adds star for significant p-val on an existing barplot.

    Args:
        ax (matplotlib.Axes): matplotlib figure (axis object) to update
        to_highlight (list of int): x values for significant bars
    """

    for i, rect in enumerate(ax.patches):
        if i in to_highlight:
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width()/2.0, height, '*', ha='center', va='bottom')


def barplot(data, output, title='', xlabel='', ylabel='', avg=None, avg_lab='', highlight_over='',
            highlight_under='', sign_all=False, sign_up_only=False, sign_down_only=False):

    """
    Plots data as a barplot and highlight significant enrichment and/or depletion. Saves the plot
    to file.

    Args:
        data (tuple of lists): categorical input data, with each tuple containing, category,
                               proportion of observed counts, enrichment or depletion, whether null
                               hyp. is rejected, and adjusted p-values.
        output (str): output file name for the figure
        title (str, optional): title for the plot
        xlabel  (str, optional): label for the x axis
        ylabel (str, optional): label for the y axis
        avg (float, optional): average over bars, plot as a dashed-line if given
        avg_lab (str, optional): label to give to the average line
        highlight_under (str, optional): color for bars under average
        highlight_over (str, optional): color for bars over average
        sign_all (bool, optional): highlight significant enrichment & depletion with stars
        sign_up_only (bool, optional): highlight significant enrichment with stars
        sign_down_only (bool, optional): highlight significant depletion with stars

    """

    if highlight_over and highlight_under:
        clrs = []
        for i in data[2]:
            if i == 'depletion':
                clrs.append(highlight_under)
            elif i == 'enrichment':
                clrs.append(highlight_over)
            else:
                clrs.append('lightgrey')

    if highlight_over:
        clrs = [highlight_over if i == 'enrichment' else 'lightgrey' for i in data[2]]
    elif highlight_under:
        clrs = [highlight_under if i == 'depletion' else 'lightgrey' for i in data[2]]
    else:
        clrs = 'lightgrey'

    bars = pd.DataFrame([(data[1][i], data[0][i]) for i in range(len(data[0]))],
                        columns=[ylabel, xlabel])

    ax = sns.barplot(x=xlabel, y=ylabel, data=bars, palette=clrs)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)

    found = False
    sys.stdout.write("### Homeologs enriched in sequence-syteny conflicts: ###\n")
    if sign_all:
        sign = [i for i, rej in enumerate(data[3]) if rej]
        for i in sign:
            sys.stdout.write(f"{data[0][i]}, pval={data[-1][i]}\n")
            found = True
    elif sign_up_only:
        sign = [i for i, rej in enumerate(data[3]) if rej and data[2][i] == 'enrichment']
        for i in sign:
            sys.stdout.write(f"{data[0][i]} significantly enriched, pval={data[-1][i]}\n")
            found = True
    elif sign_down_only:
        sign = [i for i, rej in enumerate(data[3]) if rej and data[2][i] == 'depletion']
        for i in sign:
            sys.stdout.write(f"{data[0][i]} significantly depleted, pval={data[-1][i]}\n")
            found = True
    if not found:
        sys.stdout.write("No significantly enriched homeolog.\n")


    plot_sign(ax, sign)


    if avg:
        plt.plot([-0.5, len(data[0])], [avg, avg], color='black', linestyle='--', linewidth=1)
        leg = plt.legend([avg_lab])
        # sys.stdout.write(f"{avg_lab}: {round(avg, 3)}\n")

        # as usual, seaborn does its thing...
        # Let's get the individual lines inside legend and set line width and style
        for line in leg.get_lines():
            line.set_linestyle('--')
            line.set_linewidth(1)

    if title:
        plt.title(title)

    sns.despine(trim=True)
    plt.tight_layout()
    plt.savefig(output)
    plt.close('all')


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--inc', required=True)

    PARSER.add_argument('-g', '--nb_genes', required=True)

    PARSER.add_argument('-o', '--output', type=str, required=True)

    PARSER.add_argument('--refname', type=str, required=False, default='Outgroup')

    PARSER.add_argument('-a', '--alpha', type=int, required=False, default=0.05)

    ARGS = vars(PARSER.parse_args())

    OBS = load_counts(ARGS["inc"])
    ALL = load_counts(ARGS["nb_genes"])

    for key in ALL:
        if key not in OBS:
            OBS[key] = 0

    RES = hypergeom_enrich_depl(OBS, ALL, alpha=ARGS["alpha"], multitest_adjust='fdr_bh')

    AVG_PROP = sum(OBS.values())  / sum(ALL.values()) * 100

    ALPHA = ARGS["alpha"]

    barplot(RES, ARGS["output"],
            title=f'Homeologs with increased sequence/synteny conflicts \n (* p-value < {ALPHA})',
            xlabel=f'{ARGS["refname"]} chromosomes', ylabel='Sequence/synteny conflicts (%)',
            avg_lab='Genome-wide conflicts (%)', avg=AVG_PROP, highlight_over='lightcoral',
            sign_up_only=True)
