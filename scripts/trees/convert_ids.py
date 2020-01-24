#!/usr/bin/env python

"""
    Script to convert gene IDs in the trees and alignment files to shorter IDs.
    This will allow the alignment to be converted to the phylip format so that phyml can be run
    with correct input formats (trees and ali).
    Converted output filenames are input filenames prefixed with 'tmp_'.

    Example:
        $ python -m scripts.trees.convert_ids -t gene_tree1.nh gene_tree2.nh -a ali.fa
"""

import os
import argparse

from ete3 import Tree


def convert_tree(treefile, output, d_conv=None, text=''):

    """
    Converts gene IDs in an input tree. A conversion dictionary can be given, otherwise it is
    generated.

    Args:
        treefile (file): input tree in newick format.
        output (str): name for the output file.
        d_conv (dict, optional): Conversion from old to new IDs.
        text (str, optional): Debug information

    Returns:
        dict: Conversion old to new IDs.

    """


    tree = Tree(treefile)

    if not d_conv:

        leaves = [i.name for i in tree.get_leaves()]

        #For treebest-type IDs (i.e last '_' is followed by species name):
        #generated IDs are 3 letters from species name + a unique number.
        ids = [gene.split('_')[-1][0:3]+str(nb) for nb, gene in enumerate(leaves)]
        d_conv = dict(zip(leaves, ids))

    leaves = tree.get_leaves()

    assert len(leaves) == len(d_conv), "Trees have different number of leaves {}".format(text)

    for leaf in leaves:

        assert leaf.name in d_conv, "{} present in {} but not in all trees".format(leaf.name,
                                                                                   treefile)
        leaf.name = d_conv[leaf.name]

    tree.prune([i for i in tree.get_leaves()])

    tree.write(outfile=output, format=9)

    return d_conv


def convert_ali(fastafile, output, d_conv):
    """
    Converts gene IDs in an input multiple gene alignment in fasta format.
    The conversion dictionary must be given.

    Args:
        fastafile (file): input tree in newick format.
        output (str): name for the output file.
        d_conv (dict): Conversion from old to new IDs.

    """

    nb_genes = 0

    with open(fastafile, 'r') as infile, open(output, 'w') as outfile:

        for line in infile:

            if ">" in line:

                nb_genes += 1

                leaf = line[1:-1]
                assert leaf in d_conv, "{} present in {} but not in trees".format(leaf, fastafile)
                outfile.write(line.replace(leaf, d_conv[leaf]))

            else:

                outfile.write(line)

    assert len(d_conv) == nb_genes,\
           "Trees and alignment {} have different number of genes".format(fastafile)


if __name__ == '__main__':

    ## ARGS
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--treefiles', type=str, help='Folder with constrained trees',
                        required=True, nargs='+')

    PARSER.add_argument('-a', '--alifiles', type=str, help='Folder with subalis',
                        required=True, nargs='+')

    ARGS = vars(PARSER.parse_args())

    for i, treef in enumerate(ARGS['treefiles']):

        directory, filename = os.path.split(treef)
        if directory:
            outfilename = directory + '/tmp_' + filename

        else:
            outfilename = 'tmp_' + filename

        if i == 0:

            conversion = convert_tree(treef, outfilename)

        else:

            convert_tree(treef, outfilename, conversion, ARGS['treefiles'])

    for ali in ARGS['alifiles']:

        directory, filename = os.path.split(ali)
        if directory:
            outfilename = directory + '/tmp_' + filename

        else:
            outfilename = 'tmp_' + filename
        convert_ali(ali, outfilename, conversion)
