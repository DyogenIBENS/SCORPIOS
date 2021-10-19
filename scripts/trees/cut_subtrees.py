#!/usr/bin/env python

"""
    From a synteny-derived constrained tree topology, extract genes together in an orthogroup and
    their sequence alignment, for treebest phyml independent resolution of each orthogroup.

    Example::

        $ python -m scripts.trees.cut_subtrees -t ctree.nh -a ali.fa -og outgr_gene_name
        -oa outali -ot outtree
"""

import os
import argparse
import operator

from ete3 import Tree

from .utilities import get_subali, write_fasta

def get_orthogroups_genes(ctree, outgr_gene_name):

    """
    Finds the two polytomies in the constrained tree topology.

    Args:
        ctree (str): input tree file in newick format.
        outgr_gene_name (str): gene name of the outgroup gene.

    Returns:
        dict: the 1 or 2 polytomy node(s) and their corresponding size.
        str: full outgroup gene name (with species tag)
    """

    ctree = Tree(ctree)
    orthogroups = {}
    outgr = ''

    for leaf in ctree.get_leaves():

        if outgr_gene_name != '_'.join(leaf.name.split('_')[:-1]):

            parent_node = leaf.up

            if parent_node not in orthogroups:

                orthogroups[parent_node] = len(parent_node.get_leaves())
        else:
            outgr = leaf.name

        if len(orthogroups) == 2:
            break

    return orthogroups, outgr


def write_resolved_tree(orthog_tree, outgr_gene_name, out):

    """
    Writes solution trees for orthogroup with only 2 genes.

    Args:
        orthogroup tree (ete3.Treeode) : Node with the 2 descendants of the orthogroup.
        outgr_gene_name (str): full outgroup gene name (with species tag).
        outfile (str): filename to write the tree.
    """

    new_tree = Tree()

    new_tree.add_child(orthog_tree)
    new_tree.add_child(name=outgr_gene_name)

    new_tree.prune([i for i in new_tree.get_leaves()])

    new_tree.write(outfile=out, format=1)


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-c', '--ctree', help='Newick constrained tree to resolve', required=True)

    PARSER.add_argument('-a', '--ali', help='Corresponding sub-alignment file (fasta)',
                        required=True)

    PARSER.add_argument('-og', '--outgr_gene', help='Name of the outgroup gene', required=True)

    PARSER.add_argument('-ot', '--outrees', help='Path to write resolved subtrees',
                        required=True)

    PARSER.add_argument('-oa', '--outalis', help='Path to write resolved subalis',
                        required=True)


    ARGS = vars(PARSER.parse_args())


    ORTHOGROUPS, OUTGR_GENE = get_orthogroups_genes(ARGS["ctree"], ARGS["outgr_gene"])

    #Load alignment
    with open(ARGS["ali"], 'r') as infile:
        ALI = infile.read().strip()

    #Sort orthogroups by size
    ORTHOGROUPS = sorted(ORTHOGROUPS.items(), key=operator.itemgetter(1), reverse=True)

    for k, (tree, size) in enumerate(ORTHOGROUPS):

        if size > 2:

            #write subalignment to compute a gene tree with it
            seq = get_subali(ALI, [i.name for i in tree.get_leaves()] + [OUTGR_GENE])
            write_fasta(seq, ARGS["outalis"]+'_'+str(k+1)+'.fa')

        else:
            #if 2 leaves + outgroup gene, the topology is already resolved
            outfile = ARGS["outrees"]+'_'+str(k+1)+'.nh'
            write_resolved_tree(tree, OUTGR_GENE, outfile)

            #safeguard to remove any potential artefacts (very unlikely) from previous run
            if os.path.isfile(os.path.isfile(ARGS["outalis"]+'_'+str(k+1)+'.fa')):
                os.remove(ARGS["outalis"]+'_'+str(k+1)+'.fa')

    #safeguard to remove potential artefacts (very unlikely) from previous run
    if len(ORTHOGROUPS) == 1:
        if os.path.isfile(ARGS["outalis"]+'_2.nh'):
            os.remove(ARGS["outalis"]+'_2.fa')
