#!/usr/bin/env python

"""
    This script uses synteny conservation patterns to predict orthologous gene pairs in 2
    wgd-duplicated species.

    Example:
        $ python -m scripts.synteny.pairwise_orthology_synteny -i OrthoTable.txt
                                                               -p Oryzias.latipes_Danio.rerio
                                                               -chr LG1 -ortho TreesOrthologies/
                                                               [-o out] [-w 15] [-cutoff 0]
                                                               [-filter None]
"""


import argparse
import os

from . import syntenycompare as synt
from . import utilities as ut
from . import filter_regions as freg


def load_tree_orthologies(orthology_file, rev=False):

    """
    Loads orthologies from a tabulated-separated orthology files giving pre-computed orthologous
    gene pairs in species1 and species2, based on molecular sequence evolution.

    Args:
        orthology_file (str): name of the input orthology file
        rev (bool, optional): should species in column 1 and 2 be inverted
                               (i.e use sp2 genes as dict keys)

    Returns:
        dict: For each gene in sp1 (keys), a list of orthologous genes in sp2 (values), resp. sp2
        and sp1 if rev is True.
    """

    orthos = {}

    ind1, ind2 = 0, 1

    #if rev is used column 2 (ind2) is key column 1 (ind1) is value
    if rev:
        ind1, ind2 = 1, 0

    #browse file and store orthologies
    with open(orthology_file, 'r') as infile:

        for i, line in enumerate(infile):

            line = line.strip().split('\t')

            #skip the header
            if i > 0:

                gene1 = line[ind1].strip()

                if gene1 not in orthos:
                    orthos[gene1] = []

                gene2 = line[ind2]

                if gene2 not in orthos[gene1]:
                    orthos[gene1].append(gene2)

    return orthos


def synteny_orthology_prediction(orthotable, sp1, sp2, chrom, tree_orthos, res_orthologies,
                                 win_size=15, cutoff=0, regions=None):

    """
    Compares synteny similarity of duplicated segments stored in the Orthology Table, for `sp1`
    and `sp2, using a sliding window on chromosomes `chrom` of the outgroup. Gene pairs in similar
    syntenic context are predicted orthologs.

    Args:
        orthotable (str): Name of the file with the Orthology Table
        sp1, sp2 (str): Name of compared duplicated species
        chrom (str): Name of the outgroup chromosome
        tree_orthos (dict): Orthologous gene pairs in sp1 and sp2, defined from molecular evolution
        res_orthologies (dict): dict to store results
        win_size (int, optional): Size of the sliding window to browse the orthology table
        cutoff (int, optional): cutoff on synteny similarity delta scores to predict orthology
        regions (list, optional): List of regions on the outgroup chromosome to restrict the
                                   analysis on

    Returns:
        dict: Synteny-predicted orthologous gene pairs

    """

    #load orthotable entries for sp1 and sp2
    table_entries_sp1 = ut.complete_load_orthotable(orthotable, chrom, sp1)
    table_entries_sp2 = ut.complete_load_orthotable(orthotable, chrom, sp2)

    #if no filtered regions are defined browse the whole orthology table
    if not regions:
        regions = [(0, len(table_entries_sp1))]

    #browse all windows on this outgroup chromosome in authorized regions
    for reg in regions:

        #when we load the regions the end can go out of bounds
        if reg[1] > len(table_entries_sp1):
            reg = (reg[0], len(table_entries_sp1))

        #check that regions is at least as long as windowSize
        if (reg[1] - reg[0]) + 1 >= win_size:

            #sliding window of size win_size
            for i in range(reg[0], reg[1] - win_size + 1):

                #get orthotable entry for this window and initialise a DupSegments object
                dup_seg_sp1 = synt.to_dup_segments(table_entries_sp1[i:i+win_size])
                dup_seg_sp2 = synt.to_dup_segments(table_entries_sp2[i:i+win_size])

                #if we have at least one segment in each species
                if dup_seg_sp1.matrix.shape[1] >= 1 and dup_seg_sp2.matrix.shape[1] >= 1:

                    #thread ancestrally duplicated regions
                    best, s_max = find_best_threading(dup_seg_sp1, dup_seg_sp2, tree_orthos)

                    #Store gene orthologies for best threading
                    if abs(s_max) > cutoff:
                        dup_seg_sp1.update_orthologies(dup_seg_sp2, s_max, best, res_orthologies)

    return res_orthologies


def find_best_threading(dup_seg_sp1, dup_seg_sp2, tree_orthos):

    """
    For all threading possibilities for duplicated segments in sp1 and sp2, finds the most
    parsimonious scenario.

    Args:
        dup_seg_sp1, dup_seg_sp2 (DupSegments objects): duplicated segments in sp1 and sp2
        tree_orthos (dict): Orthologous gene pairs in sp1 and sp2, defined from molecular evolution

    Returns:
        best (tuple): most parsimonious threading scenario for sp1 and for sp2
        s_max (float): corresponding synteny similarity score (delta score)
    """

    #Compute all possible threading from minor segments
    all_1 = dup_seg_sp1.all_reduce_in_two_blocks()
    all_2 = dup_seg_sp2.all_reduce_in_two_blocks()

    #Find the one that maximizes synteny similarity (delta score)
    s_max = 0
    best = (all_1[0], all_2[0])

    for red1 in all_1:

        for red2 in all_2:

            pattern_score, ortho_score = dup_seg_sp1.get_score(dup_seg_sp2, tree_orthos, red1,
                                                               red2)
            #if both scores agree, the window is informative
            if (pattern_score > 0 and ortho_score > 0) or\
               (pattern_score < 0 and ortho_score < 0):

                score = (pattern_score+ortho_score)/2.0

                if abs(score) > abs(s_max):

                    best = (red1, red2)
                    s_max = score

    return best, s_max


def write_orthologies(out, all_orthologies, sp1, sp2, filter_genes=None):

    """
    Writes synteny-predicted orthologies to file.

    Args:
        out (str): name of the output file
        all_orthologies (dict): Synteny-predicted orthologous gene pairs
        sp1, sp2 (str): Name of compared duplicated species
        filter_genes (list of str, optional): Restricted list of gene families to write (restrict
                                               the orthology prediction to some families)
    """


    with open(out, 'w') as outfile:

        #browse gene families
        for gene_family in all_orthologies:

            #filter gene families if filter_genes is set
            if not filter_genes or gene_family in filter_genes:

                score, orthologs = all_orthologies[gene_family]

                score_text = str(score)

                #for each orthology group
                for pair in orthologs:


                    #orthologs in sp2
                    genes2 = (',').join([i+'_'+sp2 for i in pair[1]])

                    #if no ortholog in sp2 write place-holder
                    if not genes2:
                        genes2 = 'None_'+sp2

                    #orthologs in sp2
                    if pair[0]:

                        for gene1 in pair[0]:

                            gene1 += '_'+sp1

                            outfile.write(gene1+'\t'+genes2+'\t'+score_text+'\t'+gene_family+'\n')

                    #if no ortholog in sp2 write place-holder
                    else:

                        gene1 = 'None_'+sp1

                        outfile.write(gene1+'\t'+genes2+'\t'+score_text+'\t'+gene_family+'\n')



if __name__ == '__main__':


    ## Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,\
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required

    PARSER.add_argument('-i', '--input', type=str, help='Orthology Table', required=True)

    PARSER.add_argument('-p', '--pair', type=str, help="Pair of species to compare,\
                         separated by '_'.", required=True)

    PARSER.add_argument('-chr', '--chr_outgr', type=str, help='Outgroup chromosome to consider.',
                        required=True)

    PARSER.add_argument('-ortho', '--orthoFromTrees', type=str,
                        help='Folder with orthologies prediction from trees.',
                        required=True)

    # Optional

    PARSER.add_argument('-o', '--out', type=str, help='Result file', required=False, default="out")

    PARSER.add_argument('-w', '--windowSize', type=int, help='Size of the sliding window',
                        required=False, default=15)

    PARSER.add_argument('-cutoff', '--cutoff', type=float, help='Synteny information cutoff',
                        required=False, default=0)

    PARSER.add_argument('-filter', '--filter_for_iter',
                        help='Filter families when using iterative mode.',
                        required=False, default=None)

    PARSER.add_argument('-chr_list', '--chr_list', type=str, help='If provided the script will run\
                        for all chromosomes of the outgroup, as provided in the file passed with\
                        this argument', required=False, default='')

    ARGS = vars(PARSER.parse_args())

    (SP1, SP2) = ARGS["pair"].split('_')

    #load pre-computed orthologies
    REV = False
    ORTHO_FILE = ARGS["orthoFromTrees"]+'/ens_'+SP1+'_'+SP2+'.txt'
    if not os.path.isfile(ORTHO_FILE):
        REV = True
        ORTHO_FILE = ARGS["orthoFromTrees"]+'/ens_'+SP2+'_'+SP1+'.txt'
    TREE_ORTHOS = load_tree_orthologies(ORTHO_FILE, REV)


    CHROMS = [ARGS["chr_outgr"]]
    if ARGS["chr_list"]:

        CHROMS = ut.outgr_chromosomes(ARGS["chr_list"])

    ALL_GENES = []

    SYNTENY_ORTHOS = {}
    for CHROM in CHROMS:

        #if SCORPiOs is run in iterative mode, we restrict regions to use in the synteny analysis
        #to regions and gene families to regions with updated synteny information
        #Here `GENES` are all genes within an updated synteny context and `REGIONS` are all windows
        #they appear in
        REGIONS, GENES = None, None
        if ARGS['filter_for_iter']:

            REGIONS, GENES = freg.read_authorized_regions(ARGS['filter_for_iter'], CHROM,
                                                          ARGS['windowSize'])
            ALL_GENES += GENES

        #predict gene orthologies using synteny
        SYNTENY_ORTHOS = synteny_orthology_prediction(ARGS["input"], SP1, SP2, CHROM,
                                                      TREE_ORTHOS, SYNTENY_ORTHOS,
                                                      ARGS["windowSize"], ARGS['cutoff'], REGIONS)

    #write predictions
    write_orthologies(ARGS["out"], SYNTENY_ORTHOS, SP1, SP2, filter_genes=ALL_GENES)
