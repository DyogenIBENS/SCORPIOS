#!/usr/bin/env python

"""
    Writes a 3-columns file for gene families, giving family_id, genes in the family and its class.
    Class can be for instance LORe or AORe, tree clustering, synteny consistency etc...

    Example::

        $ python -m scripts.lorelei.write_ancgenes_treeclust TODO
"""

import os.path
import argparse
from ete3 import Tree


def write_ancgenes(clustered_genes, treedir, out_ancgenes, clusters_to_load=None):

    """
    Writes the output 3-columns file, tab-separated.

    Args:
        clustered_genes (dict): class of gene families
        treedir (str): path to the gene trees
        out_ancgenes (str): name of the output file
        clusters_to_load (list, optional): write only entries for these given family classes.

    """

    k = 0

    with open(out_ancgenes, 'w') as outfile:

        for gene in clustered_genes:

            cluster = clustered_genes[gene]

            #Load only required family classes
            if clusters_to_load is not None and cluster not in clusters_to_load:
                continue

            #try different name for the input tree given the tree directory
            treefile = treedir +  '/' + gene + '.nhx'
            if not os.path.exists(treefile):
                treefile = treedir +  '/' + gene + '.nh'

            if not os.path.exists(treefile):
                treefile = treedir +  '/C_' + gene + '.nh'

            if not os.path.exists(treefile):
                treefile = treedir + "/" + gene + "_final.nhx"

            assert os.path.exists(treefile), f"The file {treefile} does not exist"

            tree = Tree(treefile)

            leaves = {'_'.join(i.name.split('_')[:-1]) for i in tree.get_leaves()}

            if leaves == {''}:
                leaves = {i.name for i in tree.get_leaves()}

            descendants = sorted(list(leaves))

            if clusters_to_load is not None:
                cluster = str(clusters_to_load.index(cluster))

            outfile.write(gene+'\t'+ ' '.join(descendants)+'\t'+cluster+'\n')

            k += 1


def load_gene_list(input_summary, input_acc=None):

    """
    Loads a tab-delimited summary of tree classes.

    Args:
        input_summary(str): path to the two-columns tab-delimited input file, giving a family_id to
                            tree class correspondance.
                            The family_id should be the name of the corresponding tree file for
                            write_ancgenes to work properly.
        input_acc(str, optional): if input is SCORPiOs-generated sequence-synteny inconsistent trees
                                  summary, provide here the summary of accepted correction.
                                  Indeed, gene trees that were initially found to be
                                  synteny-inconsistent but were later corrected should be defined as
                                  consistent.

    Returns:
        dict: for each gene family, the corresponding gene tree class

    """

    with open(input_summary, 'r') as infile:
        genes = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in infile}

    if input_acc is not None:

        with open(input_acc, 'r') as infile:
            acc = {line.strip().split('\t')[0] for line in infile}

        for gene in genes:
            if genes[gene] == "Inconsistent" and gene in acc:
                genes[gene] = "Consistent"

    return genes


def write_summary(summary_dict, output_file):
    """
    Writes a simpler 2-columns file with family_id and family class.

    Args:
        summary_dict (dict): class of gene families
        output_file (str): name of the output file
    """
    with open(output_file, 'w') as out:
        for key in summary_dict:
            out.write('\t'.join([key, summary_dict[key]])+'\n')


if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--treesdir', help='', required=True)

    PARSER.add_argument('-c', '--clusters', help='', required=True)

    PARSER.add_argument('-a', '--accepted', help='SCORPiOs accepted corrections',
                        required=False, default=None)

    PARSER.add_argument('-o', '--outfile', help='Output file', required=False, default="out")

    PARSER.add_argument('--summary_only', help='Only write outgroup gene name + tree consistency.',
                        action="store_true")

    PARSER.add_argument('-r', '--restrict_to', required=False, default=None, nargs='*')

    ARGS = vars(PARSER.parse_args())

    CLUSTERS = load_gene_list(ARGS["clusters"], ARGS["accepted"])

    if ARGS["summary_only"]:
        write_summary(CLUSTERS, ARGS["outfile"])

    else:
        write_ancgenes(CLUSTERS, ARGS["treesdir"], ARGS["outfile"], ARGS["restrict_to"])
