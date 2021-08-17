

"""

"""
import sys

import os
import argparse
import gzip
from collections import Counter

from ete3 import Tree

from scripts.synteny.duplicated_families import tag_duplicated_species
from scripts.trees.speciestree import get_species
import scripts.trees.utilities as ut

#TODO handle multiple outgroups better (see constrained_aore_lore)...
#small_sp_set should be dup sp only and final set small_sp + outgr
#+ allo to fetch tree named with second scorpios outgr (even if only one for clust)


def check_copy_number(tree, ref_species, sp_min=3, copy_max=2, sp_min_2copies=0, copy_in_ref=None,
                      groups=None):

    ref_gene = [i.name for i in tree.get_leaves() if i.S in ref_species]
    species = [i.S for i in tree.get_leaves()]

    if not set(species).intersection(ref_species):
        return False

    if len(set(species)) <= sp_min:
        return False

    if copy_in_ref is not None and len(ref_gene) != copy_in_ref:
        return False

    ref_gene = ref_gene[0]
    count_genes = Counter(species)

    sp_2copies = [i for i in count_genes if count_genes[i] == 2]
    if len(sp_2copies) < sp_min_2copies:
        return False

    for spec in count_genes:
        if count_genes[spec] > copy_max:
            return False

    groups_2copies = []
    if groups is not None:
        for group in groups:
            for spec in group:
                if spec in sp_2copies:
                    groups_2copies.append(spec)
                    break

        if len(groups_2copies) != len(groups):
            return False

    return True

def extract_subtrees(tree, ali, target_species, ref_species, outali, outtrees,
                     restrict_species_set=None, balanced_groups=None):

    """

    """

    tree = Tree(tree)

    #find all monphyletic groups
    tag_duplicated_species(tree.get_leaves(), target_species)

    #all clades with only target species genes in the tree
    subtrees = tree.get_monophyletic(values=["Y"], target_attr="duplicated")

    for subtree in subtrees:
        subtree_copy = subtree.copy()
        if restrict_species_set is not None:

            small_set = [i for i in subtree_copy.get_leaves() if i.S in restrict_species_set]

            if len(small_set) > 3:

                subtree_copy.prune([i for i in subtree_copy.get_leaves()\
                                    if i.S in restrict_species_set], preserve_branch_length=True)
            else:

                continue


        if not check_copy_number(subtree_copy, ref_species, sp_min=3, copy_max=2, sp_min_2copies=2,
                                 groups=balanced_groups, copy_in_ref=1):
            continue


        ref_gene = [i.name for i in subtree_copy.get_leaves() if i.S in ref_species]

        ref_gene = ref_gene[0]

        count_sp = {}

        for leaf in subtree_copy.get_leaves():

            species = leaf.S
            count_sp[species] = count_sp.get(species, 0) + 1

            leaf.S = leaf.S  + '_' + str(count_sp[species])

        subtree_copy.write(outfile=outtrees+'/'+ref_gene+'.nhx', format=1, features=["S", "D"])

        leaves = [i.name for i in subtree_copy.get_leaves()]

        seq = ut.get_subali(ali, leaves)
        ut.write_fasta(seq, outali + '/' + ref_gene + '.fa')



if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--trees', help='Forest of trees in .nhx, with species,\
                         duplication/speciation nodes and ancestor species.',\
                         required=True)

    PARSER.add_argument('-a', '--alis', help='Forest of trees in .nhx, with species,\
                         duplication/speciation nodes and ancestor species.',\
                         required=True)

    PARSER.add_argument('--anc', help='Name of the ancestor of all duplicated genomes.',\
                         required=True)

    PARSER.add_argument('-s', '--speciesTree', help='Species tree (newick), with ancestor names.',\
                         required=True)

    PARSER.add_argument('-oa', '--outdir_ali', help='Output directory for subalis',
                        required=False, default="out_alis")

    PARSER.add_argument('-ot', '--outdir_trees', help='Output directory for subtrees',
                        required=False, default="out_trees")

    PARSER.add_argument('-sp', '--outgr_species', required=True, nargs='+')

    PARSER.add_argument('-set', '--small_sp_set', required=False, default=None, nargs='+')

    PARSER.add_argument('--balance_ohno', required=False, action="store_true")

    ARGS = vars(PARSER.parse_args())

    OPEN = open

    if ARGS["alis"].split('.')[-1] == 'gz':
        OPEN = gzip.open

    #Dup. species
    SPECIES = get_species(ARGS["speciesTree"], ARGS["anc"])

    SPTREE = Tree(ARGS["speciesTree"], format=1)
    ANC = SPTREE.get_common_ancestor(SPECIES + [ARGS["outgr_species"]])
    SPECIES2 = {i.name for i in ANC.get_leaves()}

    os.makedirs(ARGS["outdir_ali"], exist_ok=True)
    os.makedirs(ARGS["outdir_trees"], exist_ok=True)

    SMALL_SP_SET = ARGS.get("small_sp_set", None)

    GROUPS = None
    if ARGS["balance_ohno"]:
        ANCESTOR = ARGS["anc"]
        SP_SET = SPECIES
        if SMALL_SP_SET is not None:
            SP_SET = [i for i in SPECIES if i in SMALL_SP_SET]
            ANCESTOR = SPTREE.get_common_ancestor(SP_SET)

        GROUP1 = [i.name for i in ANCESTOR.children[0].get_leaves() if i.name in SP_SET]
        GROUP2 = [i.name for i in ANCESTOR.children[1].get_leaves() if i.name in SP_SET]
        GROUPS = [GROUP1, GROUP2]

    k = 0

    with open(ARGS["trees"], "r") as infile_t, OPEN(ARGS["alis"], "rt") as infile_a:

        for TREE, ALI in zip(ut.read_multiple_objects(infile_t),
                             ut.read_multiple_objects(infile_a)):

            if k%500 == 0 and k != 0:
                sys.stderr.write(f"Processed {k} trees\n")

            extract_subtrees(TREE, ALI, SPECIES2, ARGS["outgr_species"], ARGS["outdir_ali"],
                             ARGS["outdir_trees"], SMALL_SP_SET, GROUPS)

            k += 1
