

"""
Script that browses scorpios corrected tree forest to build constrained AORe and LORe tree
topologies.

    Example::
        $ python -m scripts.lorelei.constrained_aore_lore_topologies TODO
"""

import sys

import os
import argparse
import gzip

from itertools import compress
from collections import Counter

from ete3 import Tree

from scripts.synteny.duplicated_families import tag_duplicated_species
from scripts.trees.speciestree import get_species
import scripts.trees.utilities as ut


#FIXME: atm OUTGROUPS, if multiple, must be a monophyletic group
#FIXME: we can simplify this code --> todo when I'll add more LORe groups

def get_species_groups(speciestree, dup_anc, outgroups, restrict_sp=None, groups_by_anc=None):

    """
    Get the 2 species groups diverging at a given speciation point (`dup_anc`) + 1 group of outgroup
    species.

    Args:
        speciestree (str): filename for the newick species tree
        dup_anc (str): speciation point to consider
        outgroups (list of str): list of outgroup species to include in the outgroup group.
        restrict_sp (list of str, optional):
        groups_by_anc (str): ancestral name for the two species groups to extract, comma separated
                             (overrules the dup_anc arg).

    Returns:
        list of str: the 3 groups of species
    """

    tree = Tree(speciestree, format=1)

    if groups_by_anc is None:
        node = tree.search_nodes(name=dup_anc)[0]
        groups = [outgroups,\
                 [i.name for i in node.children[0].get_leaves()\
                  if restrict_sp is None or i.name in restrict_sp],\
                 [i.name for i in node.children[1].get_leaves()\
                  if restrict_sp is None or i.name in restrict_sp]]

    #TODO: fix this to allow more than one group
    else:
        anc1, anc2 = groups_by_anc.split(',')
        node1 = tree.search_nodes(name=anc1)[0]
        node2 = tree.search_nodes(name=anc2)[0]
        groups = [outgroups,\
         [i.name for i in node1.children[0].get_leaves()\
          if restrict_sp is None or i.name in restrict_sp],\
         [i.name for i in node2.children[1].get_leaves()\
          if restrict_sp is None or i.name in restrict_sp]]
    return groups


def get_scorpios_aore_tree(gene_list, treefile, outgroups, outgr_gene):

    """
    Loads the AORe gene tree built by SCORPiOs.

    Args:
        gene_list (dict): dict of gene_names (key) : species_names (value) to keep in the tree
        treefile (str): name of the input tree file
        outgroups (list of str): list of outgroup species to keep/add in tree
        outgr_gene (str): name of the outgroup gene

    Returns:
        ete3.Tree : the loaded tree
    """

    tree = Tree(treefile)
    tleaves = tree.get_leaves()

    #remove sp name
    for leaf in tleaves:
        leaf.name = '_'.join(leaf.name.split('_')[:-1])

    tree.prune([i for i in tleaves if i.name in gene_list])
    leaves = {i.name for i in tree.get_leaves()}
    if leaves != set(gene_list.keys()):

        diff = set(gene_list.keys()).difference(leaves)

        outgr_node = tree.get_leaves_by_name(outgr_gene)[0]
        outgr_t = Tree()
        for gened in diff:
            if gene_list[gened] in outgroups:
                outgr_t.add_child(name=gened)
            else:
                return None #TODO: print the kind of cases covered here?
        outgr_t.add_child(name=outgr_gene)
        outgr_node.add_child(outgr_t)
    tree.prune(tree.get_leaves())

    return tree

def make_tree_from_groups(subtree_leaves, species_groups, groups_are_genes=False):

    """
    Builds a gene tree from groups of species or groups of genes.

    Args:
        subtree_leaves (list of ete3.nodes): all genes to place in the tree
        species_groups (list of str): species to group together (first group is outgroup)
        groups_are_genes (bool, optional): set to True if species_groups are groups of genes

    Returns:
        ete3.Tree : resulting gene tree
        str : one outgroup gene name, to identify the tree
    """

    tree = Tree()

    outgr, group1, group2 = species_groups

    if not groups_are_genes:

        outgr = {i.name for i in subtree_leaves if i.S in outgr}
        group1 = {i.name for i in subtree_leaves if i.S in group1}
        group2 = {i.name for i in subtree_leaves if i.S in group2}

    outgr_gene = list(outgr)[0]
    if len(outgr) >= 2:
        outgr_node = tree.add_child(name='outgr_node')
        for i in outgr:
            outgr_node.add_child(name=i)

    else:
        outgr = outgr.pop()
        tree.add_child(name=outgr)


    if group1 and group2:
        next_node = tree.add_child(name="anc_3r")
        gr1 = next_node.add_child(name="gr1")
        for i in group1:
            gr1.add_child(name=i)

        gr2 = next_node.add_child(name="gr2")
        for i in group2:
            gr2.add_child(name=i)


    elif group1:
        next_node = tree.add_child(name="anc_3r")
        for i in group1:
            next_node.add_child(name=i)

    elif group2:
        next_node = tree.add_child(name="anc_3r")
        for i in group2:
            next_node.add_child(name=i)
    tree.prune(tree.get_leaves())
    return tree, outgr_gene


def check_aore_consistent_tree(subtree, outgroups, dup_sp):

    """
    Checks that groups in synteny-consistent trees look correct and use them to build the AORe tree.

    Args:
        subtree (ete3.subtree): AORe tree
        outgroups (list of str): list of outgroup species
        dup_sp (list of str): list of duplicated species

    Returns:

        tuple: a tuple containing:
            tree (ete3.Tree): resulting gene tree, None if the tree is fishy

            outgr_gene (str) : one outgroup gene name, to identify the tree, None if the tree is
            fishy
    """

    outgr_node, node_3r = subtree.get_children()

    if outgroups[0] not in {i.S for i in outgr_node.get_leaves()}:
        node_3r, outgr_node = outgr_node, node_3r

    #No dup species in outgr
    outgr = {i.name for i in outgr_node.get_leaves() if i.S in outgroups}
    if {i for i in outgr_node.get_leaves() if i.S in dup_sp}:
        return None, None

    #No outgr in dup_species
    if {i for i in node_3r.get_leaves() if i.S in outgroups}:
        return None, None

    #extract groups below the duplication node
    if hasattr(node_3r, "D") and node_3r.D == 'Y':
        group1, group2 = node_3r.get_children()
        group1 = {i.name for i in group1.get_leaves()}
        group2 = {i.name for i in group2.get_leaves()}

    else:
        group1 = {i.name for i in node_3r.get_leaves()}
        group2 = None

    #build the AORe tree
    tree, outgr_gene = make_tree_from_groups(None, [outgr, group1, group2], groups_are_genes=True)

    return tree, outgr_gene





def check_copy_number(tree, ref_species, sp_min=3, copy_max=2, sp_min_2copies=0, copy_in_ref=None,
                      groups=None):

    """
    Checks the number of gene copies in an input tree.

    Args:
        tree (ete3.Tree): input tree
        ref_species (list of str): list of outgroup species
        sp_min (int): minimal number of species in the tree
        copy_max (int): maximal number of gene copies for any species in tree
        sp_min_2copies (int): minimal number of species with 2 gene copies
        copy_in_ref (int): number of gene copies expected in ref species
        groups (list of str): groups of species, if provided,
                              sp_min_2copies has to be verified for all groups.

    Returns:
        bool : True if criteria are met, False otherwise
    """

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


def extract_subtrees(tree, ali, target_species, ref_species, treedir, outali, olore, oaore,
                     species_groups, restrict_sp=None):

    """
    For a full gene tree, extracts subtrees and builds AORe and LORe gene tree topologies for them.
    Writes aore and lore trees to file in nhx format and corresponding multiple alignement in fasta.

    Args:
        tree (str): tree file in nhx format for the considered gene family
        ali (str): alignment fasta file for the considered gene family
        target_species (list of str): duplicated+outgroup species
        ref_species (list of str): outgroup(s) species
        treedir (str): directory with SCORPiOs constrained gene tree topologies
        outali (str): output directory for the alignment
        olore (str): output directory for the lore topology (should exist)
        oaore (str): output directory for the aore topology (should exist)
        species_groups (list of str): groups of species for the LORe topology
        restrict_sp (list of str, optional): restrict the set of duplicated species to this set
    """


    tree = Tree(tree)

    #find all monophyletic groups (clades with only target species genes in the tree)
    #called duplicated for historical reason but here I fetch outgr+dup_sp
    tag_duplicated_species(tree.get_leaves(), target_species)
    subtrees = tree.get_monophyletic(values=["Y"], target_attr="duplicated")

    for subtree in subtrees:
        subtree_copy = subtree.copy()

        if restrict_sp:
            small_set = [i for i in subtree_copy.get_leaves() if i.S in restrict_sp]
            if len(small_set) > 3:

                subtree_copy.prune([i for i in subtree_copy.get_leaves() if i.S in restrict_sp],
                                   preserve_branch_length=True)
            else:
                continue

        if not check_copy_number(subtree_copy, ref_species):
            continue

        #Build contrained AORe tree topology
        gene_list = {i.name:i.S for i in subtree_copy.get_leaves()}

        file_exist = [os.path.isfile(treedir+'/C_'+gene+".nh") for gene in list(gene_list.keys())]
        file_exist = list(compress(range(len(file_exist)), file_exist))

        if len(file_exist) == 1:
            outgr_gene = list(gene_list.keys())[file_exist.pop()]
            treefile = treedir+'/C_'+outgr_gene+".nh"
            ctree_aore = get_scorpios_aore_tree(gene_list, treefile, ref_species, outgr_gene)

        elif file_exist == []:
            dup_sp = set(target_species).difference(ref_species)
            ctree_aore, outgr_gene = check_aore_consistent_tree(subtree_copy, ref_species, dup_sp)

        else:
            continue

        #Build contrained LORe tree topology
        ctree_lore, _ = make_tree_from_groups(subtree_copy.get_leaves(), species_groups)

        #check that LORe and AORe have been succesfully built and that they are different
        if ctree_aore is not None and ctree_lore is not None:
            assert {i.name for i in ctree_lore.get_leaves()} ==\
                   {i.name for i in ctree_aore.get_leaves()}, f"{ctree_aore}, {ctree_lore}"
            comp1 = ctree_aore.compare(ctree_lore)
            comp2 = ctree_lore.compare(ctree_aore)
            comp_res = max(comp1['source_edges_in_ref'], comp2['source_edges_in_ref'])

            if comp_res != 1:

                ctree_lore.write(outfile=olore+'/'+outgr_gene+'.nh', format=9, features=["D"])
                ctree_aore.write(outfile=oaore+'/'+outgr_gene+'.nh', format=9, features=["D"])

                leaves = [i.name for i in subtree_copy.get_leaves()]

                seq = ut.get_subali(ali, leaves)

                ut.write_fasta(seq, outali + '/' + outgr_gene + '.fa')


if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--trees', help='Forest of trees in .nhx, with species,'
                        'duplication/speciation nodes and ancestor species.', required=True)

    PARSER.add_argument('-a', '--alis', help='Multiple alignments in fasta.',
                        required=True)

    PARSER.add_argument('-c', '--ctree_dir', help='Directory with scorpios constrained topologies',
                        required=True)

    PARSER.add_argument('-s', '--speciesTree', help='Species tree (newick), with ancestor names.',
                        required=True)

    PARSER.add_argument('--anc', help='Name of the ancestor of duplicated species.',
                        required=True)

    PARSER.add_argument('-sp', '--outgr_species', required=True, nargs='+')

    PARSER.add_argument('-o', '--outdir_ali', help='Output directory for subalis', required=False,
                        default="out_alis")

    PARSER.add_argument('-ol', '--outdir_lore', help='Output directory for lore ctree',
                        required=False, default="out_lore_trees")

    PARSER.add_argument('-oa', '--outdir_aore', help='Output directory for aore ctree',
                        required=False, default="out_aore_trees")

    PARSER.add_argument('-set', '--small_sp_set',
                        help="Restricted list of duplicated species to consider", required=False,
                        default=None, nargs='+')

    PARSER.add_argument('-gr', '--sp_groups', help="LORe groups", required=False,
                        default=None, nargs='+')

    ARGS = vars(PARSER.parse_args())

    OPEN = open

    if ARGS["alis"].split('.')[-1] == 'gz':
        OPEN = gzip.open

    #Dup species
    SPECIES = get_species(ARGS["speciesTree"], ARGS["anc"])

    SPTREE = Tree(ARGS["speciesTree"], format=1)
    ANC = SPTREE.get_common_ancestor(SPECIES + ARGS["outgr_species"])
    SPECIES2 = {i.name for i in ANC.get_leaves()}

    if ARGS["sp_groups"] is None:
        SP_GROUPS = get_species_groups(ARGS["speciesTree"], ARGS["anc"], ARGS["outgr_species"])

    os.makedirs(ARGS["outdir_ali"], exist_ok=True)
    os.makedirs(ARGS["outdir_aore"], exist_ok=True)
    os.makedirs(ARGS["outdir_lore"], exist_ok=True)

    SMALL_SP_SET = ARGS.get("small_sp_set", None)
    if SMALL_SP_SET is None:
        SMALL_SP_SET = SPECIES + ARGS["outgr_species"]
    else:
        SMALL_SP_SET += ARGS["outgr_species"]


    k = 0

    with open(ARGS["trees"], "r") as infile_t, OPEN(ARGS["alis"], "rt") as infile_a:

        for TREE, ALI in zip(ut.read_multiple_objects(infile_t),
                             ut.read_multiple_objects(infile_a)):

            if k%500 == 0 and k != 0:
                sys.stderr.write(f"{k} trees processed...\n")

            extract_subtrees(TREE, ALI, SPECIES2, ARGS["outgr_species"], ARGS["ctree_dir"],
                             ARGS["outdir_ali"], ARGS["outdir_lore"], ARGS["outdir_aore"],
                             SP_GROUPS, SMALL_SP_SET)

            k += 1
