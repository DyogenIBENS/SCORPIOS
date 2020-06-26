#!/usr/bin/env python

"""
    Script to load orthogroups defined in the synteny analysis, transform them into a constrained
    tree and find trees that are inconsistent with the constraints. These constrained trees will be
    saved to file, along with the corresponding original subtree and sub-alignment for later
    correction purposes.

    Example:

        $ python -m scripts.trees.inconsistent_trees -i GraphsOrthogroups -t forest_v89.nhx
                                                     -a alis_v89.fa -n Lepisosteus.oculatus
                                                     -f OrthoTable [-oc ctrees] [-oa subalis]
                                                     [-ot subtrees] [-gs GraphsCutSummary]
                                                     [-s outsummary] [-wgd '']
                                                     [-fcombin out]
"""

import os
import argparse
import sys
from collections import Counter, OrderedDict
import gzip

from ete3 import Tree

from scripts.synteny import utilities as sy
from scripts.graphs import combine_outgroups as comb
from . import utilities as ut
from . import genetree as gt



def load_pred_file(input_file, outgr, d_orthotable):

    """
    Loads predicted orthogroups after community detection in graphs.
    Stores everything in a `FamilyOrthologies` object.

    Args:
        input_file (str): Input filename
        outgr (str): Name of the corresponding outgroup species
        d_orthotable (dict): Orthologytable with outgroup-duplicated species families definition

    Returns:
        dict: For each outgroup gene (key), orthogroups in a `FamilyOrthologies` object (value).
    """

    all_fam = {}

    with open(input_file, 'r') as infile:

        for line in infile:

            line = line.strip().split('\t')
            outgr_gene = line[0][:-1]
            orthogroup = line[1:]

            if outgr_gene not in all_fam:

                all_fam[outgr_gene] = FamilyOrthologies()

            all_fam[outgr_gene].update_orthologies(outgr_gene+'_'+outgr, orthogroup)

    for outgr_gene in all_fam:

        all_fam[outgr_gene].to_constrained_tree()
        all_fam[outgr_gene].genes_in_orthotable = d_orthotable[outgr_gene]

    return all_fam


class FamilyOrthologies():

    """
    FamilyOrthologies object containing the outgroup gene, the genes in each orthogroup, all genes
    in the family in the orthologytable, and the corresponding constrained gene tree topology.
    """

    def __init__(self):

        """
        Class builder, initialized with empty objects
        """

        self.outgroup_gene = ''
        self.orthogroup_a, self.orthogroup_b = [], []
        self.genes_in_orthotable = []
        self.ctree = Tree()


    def update_orthologies(self, outgroup_gene, orthogroup):

        """
        Adds genes of one orthogroup. a's and b's are arbitrary.

        Args:
            outgroup_gene (str): name of the outgroup gene
            orthogroup (list): list of names of genes in one orthogroup
        """

        self.outgroup_gene = outgroup_gene

        if self.orthogroup_a == []:
            self.orthogroup_a = orthogroup

        else:
            self.orthogroup_b = orthogroup


    def to_constrained_tree(self):

        """
        Transforms the orthogroups + outgroup into a constrained topology, represented by an ete3
        Tree object.
        """

        #put outgroup gene as outgroup
        self.ctree.add_child(name=self.outgroup_gene)

        #add 3R duplication node
        dup_3r = self.ctree.add_child(name="dup_3r")

        #if only one orthogroup
        if self.orthogroup_a and not self.orthogroup_b:

            #add genes in group orhogroup a
            for i in self.orthogroup_a:
                i = dup_3r.add_child(name=i)

        #if two orthogroups
        elif self.orthogroup_a and self.orthogroup_b:

            #add internal node for group a
            ortho_a = dup_3r.add_child(name="orthoA")

            #add genes in group a
            for i in self.orthogroup_a:
                i = ortho_a.add_child(name=i)

            #add internal node for group b
            ortho_b = dup_3r.add_child(name="orthoB")

            #add genes in group b
            for i in self.orthogroup_b:
                i = ortho_b.add_child(name=i)


    def update_constrained_tree(self, leaves_to_place, ensembl_tree):

        """
        Adds, to the constrained tree, leaves that are under the lca in the original subtree and
        were predicted to be in the family (orthotable). These can be, for instance, genes of
        lowcov species that were discarded from the synteny analysis. They will be placed in the
        same orthogroup as its closest neighbour in the original ensembl tree.

        Args:
            leaves_to_place (list): list of the name of genes to add to the ctree.
            ensembl_tree (ete3 Tree): original gene tree.

        """

        place_in_tree = {}
        genes = [i.name for i in self.ctree.get_leaves() if i.name != self.outgroup_gene]

        #find closest neighbours in original gene tree
        while leaves_to_place:
            gene_to_place = leaves_to_place.pop()
            node = ensembl_tree.get_leaves_by_name(name=gene_to_place)[0]
            place_in_tree[gene_to_place] = gt.closest_gene_in_tree(ensembl_tree, node, genes)

        # we place all genes once all neighbours are found, to not impact the position of others
        for gene_to_place in place_in_tree:

            min_g = place_in_tree[gene_to_place].name

            #add gene in the same orthogroup as the closest gene
            sis = self.ctree.get_leaves_by_name(name=min_g)[0]
            sis.add_sister(name=gene_to_place)

            #update orthogroups
            if min_g in self.orthogroup_a:
                self.orthogroup_a.append(gene_to_place)

            elif min_g in self.orthogroup_b:
                self.orthogroup_b.append(gene_to_place)



    def is_multigenic(self):

        """
        Filters multigenic subtrees, where more duplications than just the 3R duplication is
        involved. These families are often full of errors in original gene trees and difficult to
        solve.

        Returns:
            bool: Is the subtree multigenic (True) or not (False)
        """

        is_multigenic = False

        subtrees = [self.orthogroup_a, self.orthogroup_b]
        for subtree in subtrees:

            if subtree:
                species = [i.split('_')[-1] for i in subtree]
                nb_genes_in_teleost = Counter(species)
                nbg = 0
                for key in nb_genes_in_teleost:
                    nbg += nb_genes_in_teleost[key]

                #average gene per species
                mean_gene_in_teleost = nbg/float(len(nb_genes_in_teleost))

                #do not correct if more than 1.5 gene per species and more than 2 species
                if mean_gene_in_teleost > 1.5 and len(nb_genes_in_teleost) > 2:
                    is_multigenic = True
                    break

        return is_multigenic



def get_inconsistent_trees(tree, ali, outgroups, all_families, sfile, octr, otr, oal, stats=None,
                           discard_sp=None):

    """
    For a given ensembl tree, check whether synteny-derived constrained topologies are consistent
    with it. If not, the corresponding constrained trees, ensembl subtrees and ensembl
    sub-alignments will be saved to file.

    Args:
        tree (str): ensembl tree in newick format
        ali (str): ensembl ali in fasta format
        outgroups (list): list of outgroup species used in the synteny-analysis
        all_families (dict of OrthologyFamily instances): for each outgroup genes (key) an
        OrthologyFamily instance (synteny-derived orthogroups and constrained tree topology)
        cfile (str): file to write name of synteny consistent subtrees
        mfile (str): file to write name of multigenic subtrees
        stats (dict, optional): dict to count the number of consistent and inconsistent trees

    """

    whole_tree = Tree(tree, format=1)

    for leaf in whole_tree.get_leaves():

        namesp = leaf.name + '_' + leaf.S
        leaf.prev_name = leaf.name
        leaf.name = namesp

    cached_whole_tree = whole_tree.get_cached_content(store_attr=['name', 'prev_name', 'S'])

    outgr_leaves = [i for i in cached_whole_tree[whole_tree] if i[2] in outgroups]

    #for each gene of the outgroup present in tree
    for outgr_leaf in outgr_leaves:

        #if we have a corresponding constrained tree topology
        if outgr_leaf[1] in all_families:

            ctree = all_families[outgr_leaf[1]].ctree

            cached_ctree = ctree.get_cached_content(store_attr=['name'])

            #fast way to flatten list of tuple
            ctree_leaves = set(sum(cached_ctree[ctree], ()))

            lca = whole_tree.get_common_ancestor(ctree_leaves)

            leaves = cached_whole_tree[lca]

            leaves_in_fam = [i for i in leaves\
                             if i[1] in all_families[outgr_leaf[1]].genes_in_orthotable\
                             or i[0] in ctree_leaves]

            leavesnames_in_fam = {i[0] for i in leaves_in_fam}

            #keep all genes present in the family
            if len(ctree_leaves) < len(leaves):

                to_replace_inside = leavesnames_in_fam.difference(ctree_leaves)

                all_families[outgr_leaf[1]].update_constrained_tree(to_replace_inside, lca)

            if discard_sp:
                keep = [i for i in ctree_leaves if i.split('_')[-1] not in discard_sp]
                ctree.prune(keep)
                lca = lca.copy()
                lca.prune(keep)
                leavesnames_in_fam = set(keep)

            if len(leavesnames_in_fam) <= 2:
                sfile.write(outgr_leaf[1]+"\t"+"Too few genes"+'\n')
                stats['Too few genes'] = stats.get('Too few genes', 0) + 1
                continue

            comparison = ctree.compare(lca)

            #check if the constraint is present in the tree
            if comparison['source_edges_in_ref'] != 1:


                #check if family is not too multigenic, in whih case correction is difficult
                if not all_families[outgr_leaf[1]].is_multigenic():

                    #we make a copy in case more than 1 subtree is inconsistent
                    #"newick-extended" copy is iterative (based on ete3 load/write)
                    #This way we do not risk to hit recusrion limit
                    ori_tree = lca.copy("newick-extended")
                    ori_tree.prune(leavesnames_in_fam)

                    #write original subtrees
                    ori_tree.write(outfile=otr+'/'+outgr_leaf[1]+'.nh',\
                                 format=9, features=["D"])

                    #write constrained tree topology
                    ctree.write(outfile=octr+'/C_'+outgr_leaf[1]+'.nh',
                                format=9, features=["D"])


                    #write corresponding sub-alignment
                    gene_species_mapping = dict((name, sp) for namesp, name, sp  in leaves_in_fam)
                    seq = ut.get_subali(ali, gene_species_mapping, gene_species_mapping)
                    ut.write_fasta(seq, oal + '/' + outgr_leaf[1]+'.fa')
                    sfile.write(outgr_leaf[1]+"\t"+"Inconsistent"+'\n')
                    stats['Inconsistent'] = stats.get('Inconsistent', 0) + 1

                else:
                    sfile.write(outgr_leaf[1]+"\t"+"Inconsistent_multigenic"+'\n')
                    stats['Inconsistent_multigenic'] = stats.get('Inconsistent_multigenic', 0) + 1

            else:
                sfile.write(outgr_leaf[1]+"\t"+"Consistent"+'\n')
                stats['Consistent'] = stats.get('Consistent', 0) + 1


def print_out_stats(stats_dict, wgd=''):

    """
    Prints to stdout some statistics on the number of synteny consistent subtrees.

    Args:
        stats_dict (dict): a dict counting the number of consistent and inconsistent subtrees
        wgd (str, optional): the wgd for which the synteny graphs were computed
    """

    if stats_dict:

        multi = stats_dict.get('Inconsistent_multigenic', 0)
        cons = stats_dict.get('Consistent', 0)
        incons = stats_dict.get('Inconsistent', 0)
        sm = stats_dict.get('Too few genes', 0)
        tot = cons + incons
        consp = round((cons/tot)*100, 2)
        inconsp = round((incons/tot)*100, 2)

        print('\n')
        print("------------------------TREES vs SYNTENY CONSTRAINTS------------------------")
        print(" Whole-genome duplication: {}".format(wgd))
        print("\n")
        print(" {} total subtrees with predicted synteny constraints out of {}".format(tot,
                                                                                    tot+multi+sm))
        print(" ({} discarded inconsistent multigenic subtrees)".format(multi))

        print(" {} out of {} ({} %) synteny-consistent subtrees".format(cons, tot, consp))

        print(" {} out of {} ({} %) synteny-inconsistent subtrees to correct"\
              .format(incons, tot, inconsp))

        print("----------------------------------------------------------------------------")
        print("\n")


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #### Required

    PARSER.add_argument('-i', '--input', type=str, help='File with predicted orthogroups.',
                        required=True)


    PARSER.add_argument('-t', '--treesFile', help='Forest of trees in the nhx format with species,\
                         duplication/speciation nodes + duplication confidence tags.',
                        required=True)

    PARSER.add_argument('-a', '--alisFile', help='Multiple fasta containing all multiple alignment\
                        used to build the gene trees of the forest.', required=True)


    PARSER.add_argument('-f', '--orderedFamily', type=str, help='Orthology table with family\
                        definiton.', required=True)

    PARSER.add_argument('-n', '--nonDupSp', help='Non-duplicated outgroup.', required=True)


    #### Optional

    PARSER.add_argument('-oc', '--outCons', type=str, help='out folder for constrained trees',
                        required=False, default="out_ctrees")

    PARSER.add_argument('-oa', '--outAli', type=str, help='out folder for sub-alignments',
                        required=False, default="out_subalis")

    PARSER.add_argument('-ot', '--outTree', type=str, help='out folder for subtrees',
                        required=False, default="out_subtrees")

    PARSER.add_argument('-gs', '--graphs_summary', help='summary of graph cuts, used to select\
                        best prediction when using multiple outgroups', required=False,
                        default=None)

    PARSER.add_argument('-s', '--summary', help='summary of graph cuts for multiple outgroups',
                        required=False, default='trees_sum')

    PARSER.add_argument('-wgd', '--wgd_tag', type=str,
                        help='Inform on the WGD to print it along with statistics',
                        required=False, default='')

    PARSER.add_argument('-fcombin', '--file_outgroups', type=str,
                        help='File to write summary of prediction combination across outgroups',
                        required=False, default='ouy')

    PARSER.add_argument('-di', '--discard_sp', nargs='+', default=None)

    ARGS = vars(PARSER.parse_args())


    OUTGROUPS = ARGS["nonDupSp"].split(',')

    GRAPH_CUT = OrderedDict()
    ORTHOTABLE_ALL, SUMMARY = {}, {}

    ORTHOTABLES = ARGS['orderedFamily'].split(',')

    assert len(OUTGROUPS) == len(ARGS['orderedFamily'].split(',')),\
    "inconsistent ARGS for multiple outgroups"

    assert len(OUTGROUPS) == len(ARGS['input'].split(',')),\
    "inconsistent ARGS for multiple outgroups"


    sys.stderr.write("Loading synteny graphs orthogroups...\n")

    if len(OUTGROUPS) == 1:

        ALL_GRAPH_CUT = load_pred_file(ARGS["input"], ARGS["nonDupSp"],
                                       sy.load_orthotable(ARGS['orderedFamily']))

    else:

        for k, OUTGR in enumerate(OUTGROUPS):

            #Load orthogroups in a dict gene_family : orthogr a : genes et orthogr b : genes
            GRAPH_CUT[OUTGR] = load_pred_file(ARGS["input"].split(',')[k], OUTGR,
                                              sy.load_orthotable(ORTHOTABLES[k]))

        ALL_GRAPH_CUT = comb.combine_outgroups(GRAPH_CUT, ARGS["graphs_summary"],\
                                               ARGS["file_outgroups"])


    for outfolder in [ARGS["outCons"], ARGS["outTree"], ARGS["outAli"]]:
        os.makedirs(outfolder, exist_ok=True)

    if os.path.dirname(ARGS["summary"]):

        os.makedirs(os.path.dirname(ARGS["summary"]), exist_ok=True)

    sys.stderr.write("Searching gene trees for synteny-inconsistent topologies...\n")

    STATS = {}

    OPEN = open

    if ARGS["alisFile"].split('.')[-1] == 'gz':
        OPEN = gzip.open

    with open(ARGS["treesFile"], "r") as infile_t,\
         OPEN(ARGS["alisFile"], "rt") as infile_a,\
         open(ARGS["summary"], 'w') as outfile_summary:

        for TREE, ALI in zip(ut.read_multiple_objects(infile_t),
                             ut.read_multiple_objects(infile_a)):

            get_inconsistent_trees(TREE, ALI, OUTGROUPS, ALL_GRAPH_CUT, outfile_summary,\
                                   ARGS["outCons"], ARGS["outTree"], ARGS["outAli"], STATS,\
                                   ARGS["discard_sp"])

    print_out_stats(STATS, wgd=ARGS["wgd_tag"])
