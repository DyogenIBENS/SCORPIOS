#!/usr/bin/env python

"""
    Barplots and hypergeomtric tests.

    Example:

        $ python -m scripts.lore_hunter.homeologs_tree_conflicts TODO
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
        (dict) : for each category (chromosomes), as keys, the count (value) 
    """

    data = {}
    with open(input_file, 'r') as f:
        for line in f:
            if line[0] != "#":
                counts, chrom = line.strip().split()
                try:
                    ch = int(chrom)
                except ValueError:
                    ch = chrom
                data[ch] = int(counts)
    return data


def hypergeom_enrich_depl(data_obs, data_tot, alpha=0.05, multitest_adjust='fdr_bh'):

    """
    Hypergeometric tests for enrichment and depletion with multiple testing correction.
    
    Args:
        data_obs (dict): for each category observed counts
        data_tot (dict): for each category total number of objects
        alpha (float): significance level
        multitest_adjust (str): method to adjust pvalues for

    Returns:
        (tuple of lists) : categories, corresponding proportion of observed counts, enrichment or
                           depletion, whether null hyp. is rejected, and adjusted p-values.
    """

    values = []
    direction = []
    pvalues = []
    tot = sum(data_tot.values())
    tot_obs = sum(data_obs.values())
    for ch in sorted(list(data_obs.keys())):

        exp = (tot_obs/tot)*data_tot[ch]
        obs = data_obs[ch]
        if obs > exp:
            pval = ss.hypergeom.sf(obs-1, tot, tot_obs, data_tot[ch]) #enrichment
            direction.append('enrichment')
        else:
            pval = ss.hypergeom.cdf(obs+1, tot, tot_obs, data_tot[ch]) #depletion
            direction.append('depletion')

        pvalues.append(pval)

        values.append(data_obs[ch]/data_tot[ch]*100)

    mtest = multipletests(pvalues, alpha=0.05, method=multitest_adjust)

    return (sorted(list(data_obs.keys())), values, direction, list(mtest[0]), list(mtest[1]))


def plot_sign(ax, to_highlight=None):

    for i, rect in enumerate(ax.patches):
        if i in to_highlight:
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width()/2.0, height, '*', ha='center', va='bottom')


def barplot(data, output, title='', xlabel='', ylabel='', highlight_over='', avg=None, avg_lab='',
            highlight_under='', sign_all=False, sign_up_only=False, sign_down_only=False):

    """

    Args:

    Returns:

    """
    if highlight_over:
        clrs = [highlight_over if i == 'enrichment' else 'lightgrey' for i in data[2]]
    elif highlight_under:
        clrs = [highlight_under if i == 'depletion' else 'lightgrey' for i in data[2]]
    else:
        clrs = lightgrey

    bar = pd.DataFrame([(data[1][i], data[0][i]) for i in range(len(data[0]))],
                       columns = [ylabel, xlabel])

    ax = sns.barplot(x=xlabel, y=ylabel, data=bar, palette=clrs)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)

    if sign_all:
        sign = [i for i, rej in enumerate(data[3]) if rej]
    elif sign_up_only:
        sign = [i for i, rej in enumerate(data[3]) if rej and data[2][i]=='enrichment']
        for i in sign:
            sys.stdout.write(f"{data[0][i]} significantly enriched, pval={data[-1][i]}\n")
    elif sign_down_only:
        sign = [i for i, rej in enumerate(data[3]) if rej and data[2][i]=='depletion']


    plot_sign(ax, sign)


    if avg:
        plt.plot([-0.5, len(data[0])], [avg, avg], color='black', linestyle='--', linewidth=1) #draws/tot*100
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

    ARGS = vars(PARSER.parse_args())

    observed_d = load_counts(ARGS["inc"])
    all_d = load_counts(ARGS["nb_genes"])

    res = hypergeom_enrich_depl(observed_d, all_d, alpha=0.05, multitest_adjust='fdr_bh')

    avg_prop =  sum(observed_d.values())  / sum(all_d.values()) * 100

    barplot(res, ARGS["output"],
            title='Homeologs with increased sequence/synteny conflicts \n (* p-value < 0.05)',
            xlabel=f'{ARGS["refname"]} chromosomes', ylabel='Sequence/synteny conflicts (%)',
            avg_lab='Genome-wide conflicts (%)', avg=avg_prop, highlight_over='lightcoral',
            sign_up_only=True)