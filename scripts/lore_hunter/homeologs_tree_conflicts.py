import argparse

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as ss

from statsmodels.sandbox.stats.multicomp import multipletests


## Command to generate the inputs ##

#grep LG Families/OrthoTable_Salmonidae_Esox.lucius_1 | cut -f 1 | uniq -c > nb_genes_Esox_chrom

#grep LG Families/OrthoTable_Salmonidae_Esox.lucius_1 | cut -f 1,3 > Esox_genes

#grep Incons Corrections/Trees_summary_Salmonidae_2 | cut -f 1 > salmonids_inconsistent

#grep -f salmonids_inconsistent Esox_genes | cut -f 1 | uniq -c > inconsistent_trees_by_Esox_chrom

def load_counts(input_file):
    """
    Load input files
    """
    d = {}
    with open(input_file, 'r') as f:
        for line in f:
            if line[0] != "#":
                counts, chrom = line.strip().split()
                d[chrom] = int(counts)

    norm = sum(d.values())
    # d = {k:v/norm for k,v in d.items()}
    return d, norm


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--inc', required=True)

    PARSER.add_argument('-g', '--nb_genes', required=True)

    PARSER.add_argument('-o', '--output', type=str, required=True)

    PARSER.add_argument('--refname', type=str, required=False)

    ARGS = vars(PARSER.parse_args())

    expected = {}
    enrichment = {}
    pvalues = []
    inconsistencies, draws = load_counts(ARGS["inc"])
    chr_size, tot = load_counts(ARGS["nb_genes"]) 
    # print(inconsistencies, expected, draws, tot)
    print("Genome-wide", draws/tot)

    values = {}
    for c in inconsistencies:
        print(c, inconsistencies[c]/chr_size[c])
        # print(draws/tot, chr_size[c])
        try:
            ch = int(c)
        except ValueError:
            ch = c
        values[ch] = inconsistencies[c]/chr_size[c]*100
        
        exp = (draws/tot)*chr_size[c]
        obs = inconsistencies[c]
        print(c, exp, obs)
        # exp = expected[c]*tot
        # obs = inconsistencies[c]*draws
        pval = ss.hypergeom.sf(obs-1, tot, draws, chr_size[c]) #enrichment
        p2 = ss.hypergeom.cdf(obs+1, tot, draws, chr_size[c]) #depletion
        if p2 < pval:
            pval = p2  
        pvalues.append(pval)

        print('######', ch, pval)

    pval_adj = multipletests(pvalues, alpha=0.05, method='fdr_bh')

    print(values, expected)
    print(pval_adj)
    bar = [(i, j) for (j, i) in sorted(values.items())]

    clrs = ['lightcoral' if (x > (draws/tot)*100) else 'lightgrey' for x, _ in bar]

    OUTGROUP = ARGS.get("refname", "Outgroup")

    bar = pd.DataFrame(bar, columns = ['Sequence/synteny conflicts (%)', f'{OUTGROUP} chromosomes'])

    ax = sns.barplot(x=f'{OUTGROUP} chromosomes', y="Sequence/synteny conflicts (%)", data=bar, palette=clrs)
    ax.set_xticklabels(ax.get_xticklabels(),rotation=45, ha='right', fontsize=8)
    plt.plot([-0.5, len(inconsistencies)], [draws/tot*100, draws/tot*100], color='black', linestyle='--', linewidth=1)
    leg = plt.legend(["Genome-wide conflicts (%)"])

    # as usual, seaborn does its thing... 
    # Let's get the individual lines inside legend and set line width and style
    for line in leg.get_lines():
        line.set_linestyle('--')
        line.set_linewidth(1)

    plt.title("Homeologs with increased sequence/synteny conflicts")
    sns.despine(trim=True)
    plt.tight_layout()
    plt.savefig(ARGS["output"], dpi=100)
    plt.close('all')