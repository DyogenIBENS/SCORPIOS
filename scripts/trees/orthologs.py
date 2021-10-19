#!/usr/bin/env python


"""
    Script to extract orthologous genes within a gene tree forest amongst a given list of species.
    All pairwise orthologies will be stored in the output folder (one file for each species pair).

    Example::

        $ python -m scripts.trees.orthologs -t gene_trees.nhx -d Clupeocephala -s sptree.nwk
        [-o out] [-ow Salmonids] [-l lowcov_sp1,lowcov_sp2]
"""

import itertools
import os
import sys
import argparse
import errno

from ete3 import Tree

from . import utilities as ut
from . import speciestree  as spt


def is_speciation(node):

    """
    Is the node a speciation node?

    Args:
        tree (ete3.TreeNode): input node, with duplications annotated with the `D` attribute.
                              D=Y if duplication, D=N otherwise. Note that dubious nodes
                              (DD=Y or DCS=0) are considered speciation nodes.

    Returns:
        bool: True if speciation, False otherwise.

    """

    speciation = False

    if (hasattr(node, "D") and node.D == 'N'):
        speciation = True

    elif (hasattr(node, "DD") and node.DD == 'Y'):
        speciation = True

    elif (hasattr(node, "DCS") and float(node.DCS) == 0.0):
        speciation = True

    return speciation


def get_speciation_events(tree, species_pairs, sp_ortho_dict):

    """
    Extracts all orthologies relationships in a gene tree involving the species given in input,
    and adds them to the orthology dict.

    Args:
        tree (ete3.Tree): input Tree object.
        species_pairs (list): species pairs to consider.
        sp_ortho_dict (dict): dictionary to store orthologies.
    """

    #browse the tree
    for node in tree.traverse():

        #ignore leaves or artefactual single child internal nodes
        if len(node.get_children()) == 2:

            #get type of event according to duplication annotations, and continue if speciation
            #'dubious' nodes i.e duplications with confidence 0 are considered speciations
            if is_speciation(node):

                children = node.get_children()

                #get leaves separated by the speciation
                side_a_leaves = {i for i in children[0].get_leaves()}
                side_b_leaves = {i for i in children[1].get_leaves()}

                #get species at both sides of event
                inspecies = {i.S for i in side_a_leaves}
                outspecies = {i.S for i in side_b_leaves}

                #store genes separated by speciation in all species pairs
                for sp_pair in species_pairs:
                    sp1, sp2 = sp_pair
                    if (sp1 in inspecies and sp2 in outspecies) or (sp1 in outspecies and sp2 in\
                        inspecies):

                        if sp1 in inspecies:
                            genes_1 = {i.name for i in list(side_a_leaves) if i.S == sp1}
                            genes_2 = {i.name for i in list(side_b_leaves) if i.S == sp2}

                        else:
                            genes_1 = {i.name for i in list(side_b_leaves) if i.S == sp1}
                            genes_2 = {i.name for i in list(side_a_leaves) if i.S == sp2}

                        sp_ortho_dict[(sp1, sp2)] = sp_ortho_dict.get((sp1, sp2), [])

                        for (gene1, gene2) in itertools.product(genes_1, genes_2):
                            sp_ortho_dict[(sp1, sp2)].append((gene1, gene2))



if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)


    #Required
    PARSER.add_argument('-t', '--treesFile', help='Forest of trees in newhampshire format (nhx),\
                        with species, duplication/speciation nodes + duplication confidence tags.',
                        required=True)

    PARSER.add_argument('-d', '--dupSp', help='Name of the ancestor of duplicated species.',
                        required=True)

    PARSER.add_argument('-s', '--speciesTree', help='Species tree (newick), with ancestor names.',
                        required=True)

    #Optional
    PARSER.add_argument('-o', '--out', help='Result folder.', required=False, default="out")

    PARSER.add_argument('-ow', '--other_wgds', help='Ancestor(s) to exclude (Comma-delimited,\
                        exclude all species below these ancestors).', required=False, default='')

    PARSER.add_argument('-l', '--lowcov', type=str, help='Species to exclude (Comma-delimited).',
                        required=False, default='')

    ARGS = vars(PARSER.parse_args())

    #Get all pairs of duplicated species
    SPECIES = spt.get_species(ARGS["speciesTree"], ARGS["dupSp"], ARGS["other_wgds"],
                              ARGS["lowcov"])

    SPECIES = list(itertools.combinations(SPECIES, 2))

    ORTHOLOGIES = {}

    sys.stderr.write("Browsing gene trees for orthologies between duplicated species...\n")

    #browse gene tree forest
    with open(ARGS["treesFile"], "r") as infile:

        #for each gene tree
        for TREE in ut.read_multiple_objects(infile):
            TREE = Tree(TREE, format=1)

            #find all speciation nodes and corresponding orthologies
            get_speciation_events(TREE, SPECIES, ORTHOLOGIES)

        sys.stderr.write("OK\n")

    #create output folder if it does not exist
    OUT_PATH = ARGS["out"]
    try:
        os.mkdir(OUT_PATH)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise

    #Write orthologies to result files
    sys.stderr.write("Writing orthologies between duplicated species...\n")
    for (SP1, SP2) in ORTHOLOGIES:
        with open(OUT_PATH+'/ens_'+SP1+'_'+SP2+'.txt', 'w') as OUTFILE:
            for (GENE1, GENE2) in ORTHOLOGIES[(SP1, SP2)]:
                OUTFILE.write(GENE1+'\t'+GENE2+'\n')
