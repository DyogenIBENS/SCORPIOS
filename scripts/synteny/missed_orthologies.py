#!/usr/bin/env python

"""
    This script finds potential orthologs between an outgroup and duplicated species based
    on synteny, for genes without obvious orthologs in trees.

    Example::

        $ python -m scripts.synteny.missed_orthologies -i Orthotable -u UncertainGenes -c Chroms
        [-o output] [-wgd ''] [-w 0] [-f out]
"""

import collections
import argparse
import sys
import pickle

from . import utilities as syu
from . import mygenome as myg


def load_genes(genes, outgr=False):

    """
    Parses an entry in the "no phylogenetic ortholog" file and loads genes as GeneSpeciesPosition
    namedTuples.

    Args:
        genes (str): a line of the input file
        outgr (bool): Whether entry of ingroups (True) or outgroup should be parsed (False)

    Returns:

        tuple: a tuple containing:
            dict_genes (dict): for each species (key), genes in the entry (value) as a
            `GeneSpeciesPosition` namedtuple

            unplaced_genes (dict): stores genes with no gene position entry in the .bed file in a
            dict of similar structure as dict_genes
    """

    dict_genes = {}
    unplaced_genes = {}
    genes = genes.split(' ')
    for gene in genes:
        all_gene_info = gene.split('|')

        if outgr:

            name, chrom, position, _, _ = all_gene_info

            if 'Outgroup' not in dict_genes:
                dict_genes['Outgroup'] = []

            dict_genes['Outgroup'].append(syu.GeneSpeciesPosition(name, chrom,\
                                                                     myg.toint(position)))

        else:
            name, chrom, position = all_gene_info
            species = name.split('_')[-1]
            name = name.replace('_'+species, '', 1)

            if species not in dict_genes:
                dict_genes[species] = []

            if species not in unplaced_genes:
                unplaced_genes[species] = []

            if chrom:
                dict_genes[species].append(syu.GeneSpeciesPosition(name, chrom,\
                                                                      myg.toint(position)))
            else:
                unplaced_genes[species].append(syu.GeneSpeciesPosition(name, chrom,\
                                               myg.toint(position)))

    return dict_genes, unplaced_genes


def search_closest_neighbours(ingroup_genes, dup_sp, all_genefam, all_outgroup_candidates):

    """
    Extracts orthologs, in the outgroup species, of genes in the neighbourhood of genes without
    phylogenetic orthologs in species `dup sp` .

    Args:
        ingroup_genes (dict): a clade of ingroup genes without phylogenetic orthologs, as a dict,
                              giving, for each species, a list of `GeneSpeciesPosition` tuples.

        dup_sp (str): name of the considered duplicated species.

        all_genefam (nested dict): Pre-computed orthology table based on phylogenetic orthologs
                                   used to search for syntenic neighbours, represented by a nested
                                   dict, giving for each outgroup chromosome (key1) and each
                                   duplicated species (key2), a list of `GeneFamily` objects.

        all_outgroup_candidates (list) : all outgroup genes in the same tree, as a list of
                                         `GeneSpeciesPosition` tuples.

    Returns:

        tuple: a tuple containing:
            ortho_neighbours (list): a list of orthologs of ingroup genes in the outgroup,
            in a tuple (chromosome, gene index)

            skip (bool): If True, we should not use `dup_sp` to search for syntenic neighbours
            because one neighbour is orthologous to another outgroup gene in the same tree (i.e
            history of tandem duplication which will artefactually inflate the number
            of syntenic neighbours). Conservation of synteny in the case of tandem
            duplication is not a proof for orthology.
    """

    ortho_neighbours = []
    skip = False
    all_outgroup_candidates = [i.name for i in all_outgroup_candidates]
    for gene in ingroup_genes[dup_sp]:
        for outgr_chrom in all_genefam:
            for family in all_genefam[outgr_chrom][dup_sp]:
                if gene.chromosome in family.involved_chromosomes:

                    neighbour_genes = family.all_duplicate_genes
                    for nei in neighbour_genes:
                        if nei.chromosome == gene.chromosome:
                            dist = abs(nei.index - gene.index)
                            if dist <= 15 and dist != 0:
                                ortho_neighbours.append((outgr_chrom, family.outgr_position))

                                if family.outgr_genename in all_outgroup_candidates:
                                    skip = True #if neighbour of gar gene also in tree: tandem dup

    return ortho_neighbours, skip


def neighbour_outgr_ortholog(ortho_neighbours, all_outgroup_candidates):

    """
    Searches for syntenic neighbours between ingroup and outgroup genes.
    Gene neighbouring ingroup genes have their orthologs in the outgroup stored in
    `ortho_neighbours`. This function searches if `ortho_neighbours` are in the neighbourhood of an
    outgroup gene `all_outgroup_candidates` (genes in the same tree as ingroup genes).

    Args:
        ortho_neighbours (list): list of orthologs of neighbours of ingroup genes, as tuples
                                (chromosome, index)
        all_outgroup_candidates (list) : all outgroup genes in the same tree, as a list of
                                         `GeneSpeciesPosition` tuples.

    Returns:
        list: list of outgroup genes in the same tree with at least one syntenic neighbour, with
        repetitions. The number of repetitions indicates the number of syntenic neighbours.
        For instance, [gene_a, gene_a, gene_b, gene_a, gene_a] indicates that gene a has
        four syntenic neighbours with ingroup genes and gene_b one.

    """
    outgroup_genes = []
    for outgr_gene in all_outgroup_candidates:
        for (outgr_chr, outgr_pos) in ortho_neighbours:
            if outgr_gene.chromosome == outgr_chr:
                if abs(outgr_gene.index - outgr_pos) <= 15: #15 genes on both sides, fixed value
                    outgroup_genes.append(outgr_gene)
    return outgroup_genes


def find_synteny_orthologs(input_file, optimize=False, threshold=2.0, opt_fam=None):

    """
    Browses ingroup genes without phylogenetic orthologs in the outgroup and attempts to find
    synteny-supported orthologs

    Args:
        input_file (str): name of the input file storing genes without orthologs in ingroups

        optimize (bool, optional): option to use if the script is called to optimize the threshold

        threshold (float, optional): synteny support threshold

        opt_fam (list, optional): if defined, restricts fmailies to use for optimization to the
                                  ones in this list

    Returns:
        dict: identified synteny-supported orthologies, stored in nested dict with, for each
        outgroup gene with newly identified ortholog(s) (`GeneSpeciesPosition` tuple, key1)
        and for each duplicated species with such ortholog(s) (str, key2), orthologous gene
        as `GeneSpeciesPosition` tuple(s).
    """

    with open(input_file, 'r') as infile:

        sys.stderr.write("Searching for synteny-supported orthologies\n")

        res = {}

        store_scores = []

        fam_store = []

        for line in infile:

            line = line.strip().split('\t')

            if len(line) == 3:

                ingroup_genes_raw, outgroup_genes, outgroup_inparalogs = line
                outgroup_inparalogs = [tuple(sorted(i.split('|')))\
                                       for i in outgroup_inparalogs.split()]

            else:

                ingroup_genes_raw, outgroup_genes = line
                outgroup_inparalogs = []

            if opt_fam and ingroup_genes_raw not in opt_fam:
                continue

            #load clade of duplicated species genes
            ingroup_genes, unplaced_genes = load_genes(ingroup_genes_raw)

            #load all outgroup homologous genes
            outgroup_candidates = load_genes(outgroup_genes, outgr=True)[0]['Outgroup']

            synteny_outgroup_genes = []

            # for each gene in dup spec, search orthologies between gar and its neighbours
            for spec in ingroup_genes:

                #extract neighbours
                ortho_neighbours, skip = search_closest_neighbours(ingroup_genes, spec, ALL_ORTHOS,
                                                                   outgroup_candidates)

                #if no small tandem duplication in considered dup species
                #(would artificially inflate number of neighbours)
                if not skip:

                    #Are there syntenic neighbours in the 30-window centred on one duplicated gene
                    #and the 30-window centred on a gar gene present in the tree?
                    #If yes then there is synteny evidence for orthology
                    synteny_outgroup_genes += neighbour_outgr_ortholog(ortho_neighbours,
                                                                       outgroup_candidates)

            #Store outgroup gene(s) if synteny evidence found
            #Choose the one with most syntenic neighbours, if several possible.
            if synteny_outgroup_genes:

                old_tandem = False

                data = collections.Counter(synteny_outgroup_genes)

                #sort counter to find outgroup gene with the highest synteny evidence
                sorted_data = sorted(data.items(), key=(lambda x: x[1]))
                best, max_value = sorted_data[-1]

                if len(sorted_data) > 1:
                    second_best, second_max = sorted_data[-2]


                    #if the two top outgroup genes have same synteny scores --> tandem duplicate
                    #if this two genes are not recent duplications
                    #--> we cannot chose which is the most likely ortholog based on synteny
                    if max_value == second_max and\
                       tuple(sorted((best.name, second_best.name))) not in outgroup_inparalogs:
                        old_tandem = True


                #Cut-off for synteny: orthology evidence if on average per species more than 2
                #neighbour genes supports orthololgy
                #we could think of normalizing by the number of genes rather than on species:
                #the up would be to avoid over-aggregation of big clades in case of very wrong
                #trees
                #the down would be that poorer assembly would likely make the score drop
                #artificially (more gene splits --> higher number of genes in small scaffolds)
                if (max_value / float(len(ingroup_genes)) >= threshold and not optimize)\
                   and not old_tandem:


                    #store orthology
                    res[best] = res.get(best, {})

                    for spec in ingroup_genes:

                        res[best][spec] = res[best].get(spec, [])
                        res[best][spec] += ingroup_genes[spec]
                        res[best][spec] += unplaced_genes[spec]

                elif optimize and not old_tandem:
                    store_scores.append(max_value / float(len(ingroup_genes)))
                    fam_store.append(ingroup_genes_raw)

    return res, store_scores, fam_store

def print_out_stats(stats_dict, wgd='', file_fam_nograph='out_nog'):

    """
    Prints to stdout some statistics on the families in the final Orthology Table.

    Args:
        stats_dict (dict): a dict counting number of families and genes in the families

        wgd (str, optional): the wgd for which the Orhtology Table was built

        file_fam_nograph (str, optional): file to write families that can't result in a graph (won't
                          be in a large enough window or has too few genes)

    """

    empty_table = True

    if stats_dict:

        fam = stats_dict.get('families', 0)
        fam1 = stats_dict.get("families_1g", [])
        fam2 = stats_dict.get("small_chr", [])

        famg = fam - len(fam1) - len(fam2)

        if fam:

            empty_table = False

            print('\n')
            print("---------------------------FINAL ORTHOLOGY TABLE----------------------------")
            print(" Whole-genome duplication: {}".format(wgd))
            print("\n")
            print(" {} total families in the final orthology table".format(fam))
            for spec in stats_dict:
                if spec not in ['families_1g', 'families', 'small_chr']:
                    genes = stats_dict[spec]
                    print(" {} {} genes in the final orthology table".format(genes, spec))
            print("\n")
            print(" For {} families, a synteny orthology graph can potentially be built"\
                    .format(famg))
            if fam1 or fam2:
                print(" ({} with too few genes & {} in a too short window)"\
                      .format(len(fam1), len(fam2)))
                with open(file_fam_nograph, 'w') as outfile:
                    for family in fam1+fam2:
                        outfile.write(family+'\n')
            print("----------------------------------------------------------------------------")
            print("\n")

    if empty_table:
        print('\n')
        print("---------------------------FINAL ORTHOLOGY TABLE----------------------------")
        print("Whole-genome duplication: {}".format(wgd))
        print("\n")
        print("0 total families in the final orthology table")
        print("----------------------------------------------------------------------------")
        print("\n")


if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,\
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required

    PARSER.add_argument('-i', '--input', help='Orthology table with the outgroup.', required=True)

    PARSER.add_argument('-u', '--uncertain_families', help='File with genes without ortholog in\
                        their outgroup and potential orthologs.', required=True)

    PARSER.add_argument('-c', '--chr', help='File with a list of chromosomes in the outroup.',
                        required=True)

    # Optional
    PARSER.add_argument('-o', '--output', help='Output file.', required=False, default="out")

    PARSER.add_argument('-wgd', '--wgd_tag', type=str,
                        help='Inform on the corrected WGD to print it out along with statistics',
                        required=False, default='')

    PARSER.add_argument('-w', '--windowsize', type=int,
                        help='Inform on the windowsize that will be used in the synteny analysis,\
                              to print it out statistics on the number of families that will be\
                              considered in windows',
                        required=False, default=0)

    PARSER.add_argument('-f', '--fam_nograph', type=str,
                        help='File to write families that cannot result in a graph, so that it can\
                              be read by the script that filter un-updated regions. This allows\
                              to print out statistics on families in updated regions',
                        required=False, default='out_nog')

    PARSER.add_argument('-s', '--support', type=float,
                        help='Optimal number of average syntenic orthologs to recover orthologs',
                        required=False, default=2.0)


    PARSER.add_argument('-opt', '--optimize',
                        help='Run the script in optimize mode: the distribution of the number of\
                              syntenic orthologs between ingroup genes and the reference outgroup\
                              will be returned. No orthology table will be created.',
                        required=False, dest='optimize', action='store_true')

    PARSER.add_argument('-opt_fam', '--optimized_fam',
                        help='Write name of families used in the optimization',
                        type=str, default="")

    PARSER.add_argument('-u_opt_fam', '--use_optimized_fam',
                        help='Use families in input to optimize',
                        type=str, default="")

    PARSER.set_defaults(optimize=False)
    PARSER.set_defaults(optimized_fam=False)
    ARGUMENTS = vars(PARSER.parse_args())


    #Extract header with species list
    with open(ARGUMENTS["input"], 'r') as myinfile:
        for LINE in myinfile:
            OUTGR = LINE.strip().split('\t')[0]
            SP_LIST = LINE.strip().split('\t')[4:]
            break

    sys.stderr.write("Loading phylogenetic orthologies\n")

    #Load orthology table with phylogenetic orthologies
    OUTGR_CHROM = syu.outgr_chromosomes(ARGUMENTS["chr"])
    ALL_ORTHOS = collections.defaultdict(dict)
    for CHROM in OUTGR_CHROM:
        for SPECIES in SP_LIST:
            ALL_ORTHOS[CHROM][SPECIES] = syu.complete_load_orthotable(ARGUMENTS["input"], CHROM,\
                                                                      SPECIES,
                                                                      load_no_position_genes=True)
    FAM = None
    if ARGUMENTS["use_optimized_fam"]:
        with open(ARGUMENTS["use_optimized_fam"]) as myinfile:
            FAM = [line.strip() for line in myinfile]

    #Searches for synteny supported orthologies
    RES, SCORES, FAM = find_synteny_orthologs(ARGUMENTS["uncertain_families"],
                                              ARGUMENTS["optimize"],
                                              ARGUMENTS["support"],
                                              FAM)

    if ARGUMENTS["optimize"]:


        OUT = ARGUMENTS["uncertain_families"] + "_scores.pkl"
        with open(OUT, 'wb') as OUTFILE:
            pickle.dump(SCORES, OUTFILE)

        if ARGUMENTS["optimized_fam"]:
            with open(ARGUMENTS["optimized_fam"], 'w') as OUTFILE:
                OUTFILE.write('\n'.join(FAM))

    else:

        #Update the orthology table
        sys.stderr.write("Insertion of synteny-supported orthologs in the orthology table\n")
        syu.update_orthologytable(ALL_ORTHOS, RES, SP_LIST)

        #Write the new orthology table
        sys.stderr.write("Writing updated Orthology Table\n")
        STATS = syu.write_updated_orthotable(ALL_ORTHOS, OUTGR, SP_LIST, OUTGR_CHROM,
                                             ARGUMENTS["output"], wsize=ARGUMENTS['windowsize'])

        if ARGUMENTS["wgd_tag"]:
            WGD, OUTGR = ARGUMENTS["wgd_tag"].split(',')
            WGD_TAG = WGD + ' (outgroup '+OUTGR+')'
            print_out_stats(STATS, wgd=WGD_TAG, file_fam_nograph=ARGUMENTS['fam_nograph'])
