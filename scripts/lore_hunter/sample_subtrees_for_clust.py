

"""

"""
import sys

import os
import argparse
import gzip
from ete3 import Tree

from collections import Counter

from scripts.synteny.duplicated_families import tag_duplicated_species
from scripts.trees.speciestree import get_species
import scripts.trees.utilities as ut


def extract_subtrees(tree, ali, target_species, ref_species, outali, outtrees,
                     restrict_species_set=None, balanced_groups=None):

    """

    """

    tree = Tree(tree)

    #find all monphyletic groups
    tag_duplicated_species(tree.get_leaves(), target_species)

    #all clades with only target species genes in the tree
    subtrees = tree.get_monophyletic(values=["Y"], target_attr="duplicated")

    for subtree in subtrees: #ENSPKIG00000000415
        subtree_copy = subtree.copy()
        if restrict_species_set is not None:

            small_set = [i for i in subtree_copy.get_leaves() if i.S in restrict_species_set]

            if len(small_set) > 3:
            
                subtree_copy.prune([i for i in subtree_copy.get_leaves()\
                                    if i.S in restrict_species_set], preserve_branch_length=True)
            else:

                continue

        ref_gene = [i.name for i in subtree_copy.get_leaves() if i.S == ref_species]

        species = [i.S for i in subtree_copy.get_leaves()]

        if ref_species not in species:
            continue

        count_genes = Counter(species)
        sp_2copies = [sp for sp in count_genes if count_genes[sp]==2]
        if len(sp_2copies) < 2:
            continue

        if balanced_groups is not None:
            ok = False
            for i in balanced_groups[0]:
                if i in sp_2copies:
                    ok = True
                    break

            if not ok:
                continue

            ok = False
            for i in balanced_groups[1]:
                if i in sp_2copies:
                    ok = True
                    break

            if not ok:
                continue

        ok = True
        # print(subtree)
        # print(sp_2copies)

        for sp in count_genes:

            if sp == ref_species and count_genes[sp] != 1:

                ok = False
                break

            elif count_genes[sp] > 2:

                ok = False
                break


        if ok and len(set(species)) > 2:

            ref_gene = ref_gene[0]

            count_sp = {}

            for leaf in subtree_copy.get_leaves():

                sp = leaf.S
                count_sp[sp] = count_sp.get(sp, 0) + 1

                leaf.S = leaf.S  + '_' + str(count_sp[sp])

            subtree_copy.write(outfile=outtrees+'/'+ref_gene+'.nhx',
                               format=1, features=["S", "D"])

            leaves = [i.name for i in subtree_copy.get_leaves()]

            d = []
            for leaf in leaves:
                d.append(leaf) 

            seq = ut.get_subali(ali, d)
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

    PARSER.add_argument('-oa', '--outdir_ali', help='Output directory for subalis', required=False,
                        default="out_alis")

    PARSER.add_argument('-ot', '--outdir_trees', help='Output directory for subtrees', required=False,
                        default="out_trees")

    PARSER.add_argument('-sp', '--outgr_species', required=False, default=None)

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
        ancestor = ARGS["anc"]
        SP_SET = SPECIES
        if SMALL_SP_SET is not None:
            SP_SET = [i for i in SPECIES if i in SMALL_SP_SET]
            ancestor = SPTREE.get_common_ancestor(SP_SET)

        GROUP1 = [i.name for i in ancestor.children[0].get_leaves() if i.name in SP_SET]
        GROUP2 = [i.name for i in ancestor.children[1].get_leaves() if i.name in SP_SET]
        GROUPS = [GROUP1, GROUP2]

    k = 0

    with open(ARGS["trees"], "r") as infile_t, OPEN(ARGS["alis"], "rt") as infile_a:

        for TREE, ALI in zip(ut.read_multiple_objects(infile_t),
                             ut.read_multiple_objects(infile_a)):

            if k%500 == 0 and k != 0:
                sys.stderr.write(f"Processed {k} trees\n")

            extract_subtrees(TREE, ALI, SPECIES2, ARGS["outgr_species"], ARGS["outdir_ali"],
                             ARGS["outdir_trees"], SMALL_SP_SET, GROUPS)

            k+=1
