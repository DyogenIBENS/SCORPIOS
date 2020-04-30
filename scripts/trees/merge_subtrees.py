#!/usr/bin/env python

"""
    Script to merge together independently resolved orthogroups of the same family into a single
    tree.

    Example:
        $ python -m scripts.trees.merge_subtrees -t orthogroup_tree1.nh orthogroup_tree2.nh
                                                 -outgr gene_name [-o out]
"""


import argparse

from ete3 import Tree


def remove_outgroup(tree, outgr):

    """
    Loads a subtree and removes the outgroup gene.

    Args:
        tree (ete3.Tree): Input trree
        outgr (str): Outgroup gene name

    """
    tree = Tree(tree)
    leaves = [i.name for i in tree.get_leaves()]

    outgr_gene = [i for i in leaves if outgr == '_'.join(i.split('_')[:-1])][0]
    tree.set_outgroup(tree&outgr_gene)

    tree.prune([i for i in leaves if i != outgr_gene])
    return tree, outgr_gene


def merge_trees_and_write(trees, outgr, outfile, keep_br=False):

    """
    Merges two subtrees independently resolved into a single tree and adds the outgroup gene.
    Writes the result to file.

    Args:
        trees (list of ete3.Tree): Tree(s) to merge
        outgr (str): Outgroup gene name
        outfile (str): Output filename
    """

    merged_tree = Tree()

    for tree in trees:
        merged_tree.add_child(tree)

    #merge the two and place outgroup correctly
    merged_final = Tree()
    merged_final.add_child(merged_tree)
    merged_final.add_child(name=outgr)
    merged_final.prune([i for i in merged_final.get_leaves()])

    if keep_br:
        merged_final.write(outfile=outfile)

    else:
        merged_final.write(outfile=outfile, format=9)


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-outgr', '--outgroup', help='Outgroup gene', required=True)
    PARSER.add_argument('-t', '--trees', help='Resolved binary subtrees.', required=True,
                        nargs='+')
    PARSER.add_argument('-o', '--outfile', required=False, default='out')
    PARSER.add_argument('-br', '--brlength', action='store_true')
    ARGS = vars(PARSER.parse_args())

    OUTGR = ARGS["outgroup"]

    TREES = []
    for i, subtree in enumerate(ARGS['trees']):
        if i < 2:
            subtree, outgr_genename = remove_outgroup(subtree, OUTGR)
            TREES.append(subtree)
        else:
            raise ValueError("Only a maximum of two trees are expected")

    merge_trees_and_write(TREES, outgr_genename, ARGS['outfile'])
