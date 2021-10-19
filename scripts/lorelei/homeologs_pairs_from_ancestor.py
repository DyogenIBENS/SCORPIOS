"""
Loads gene tree classes and group them by ancestrally duplicated chromsome pairs.

An ancestral karyotype or a non-duplicated outgroup can be used as a proxy to the ancestral
pre-duplication genome.

    Example::

        $ python -m scripts.lorelei.homeologs_pairs_from_ancestor.py TODO
"""

import sys
from collections import Counter
import argparse

from scripts.synteny.mygenome import Genome


def load_pm(input_file):

    """
    Loads predicted homeolog names for teleosts gene families.

    Args:
        input_file (str): path to the input file (output from the paralogy_map pipeline)

    Returns:
        dict: for each gene family (key, a set of teleost genes) its corresponding homoelog
        chromosome (value).
    """

    fam_homeo = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            _, genes, homeo = line.strip().split('\t')
            genes = frozenset(genes.split())
            if homeo != "?":
                fam_homeo[genes] = int(homeo[:-1])
    return fam_homeo

def load_summary(input_file, accepted):

    """
    Loads SCORPiOs summary of synteny-sequence trees inconsistencies.

    Args:
        input_file (str): path to the input file
        accepted (str): path to the file with accepted correction (inconsistent trees which have
                        been corrected are now consistent)

    Returns:
        dict: for each outgroup gene family identifier (key, str) wheter trees and synteny
        predictions are consistent or inconsistent (value, str).
    """

    d_summary = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            fam_id, consistency = line.strip().split('\t')
            consistency = consistency.lower().split('_')[0]

            if consistency in ["consistent", "inconsistent"]:
                if consistency == "inconsistent" and fam_id in accepted:
                    consistency = "consistent"
                d_summary[fam_id] = consistency

    return d_summary


def load_outgr_fam(input_file, ctrees=None):

    """
    Loads SCORPiOs teleost families file.

    Args:
        input_file (str): path to the input file
        ctrees (set, optional): list of families to load, by default everything is loaded.

    Returns:
        dict: for each family, identified by the outgroup gene (key, str), teleost genes
        (value, set)
    """

    d_fam = {}
    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):
            if i != 0: #skip header
                line = line.strip().split('\t')
                outgr_gene = line[2]
                if ctrees is None or outgr_gene in ctrees:
                    genes = {j.split('|')[0] for gene in line[3:] for j in gene.split('/')}
                    d_fam[outgr_gene] = genes

    return d_fam

def load_acc(input_file):

    """
    Loads accepted correction.

    Args:
        input_file (str): path to the input file.

    Returns:
        set: all family ids (outgroup gene name) for which correction was accepted
    """

    with open(input_file, 'r') as infile:
        res = {line.strip().split('\t')[0] for line in infile}
    return res


def outgroup_genes_to_homeologs(fam_outgr, fam_homeo):

    """
    Combines teleost genes in SCORPiOs families and paralogy map result to assign SCORPiOs families
    to homeologs.

    Args:
        fam_outgr (dict): for each gene family (name of the outgroup gene, str, key) the set of
                          genes (value)
        fam_homeo (dict): for each gene family (set of genes, set, key) its homeolog (value)

    Returns:
        dict: for each gene family (name of the outgroup gene, str, key), its homeolog (value)
    """

    d_homeo = {}
    i = 0
    for outgr, fam1 in fam_outgr.items():

        if i % 500 == 0 and i != 0:
            sys.stderr.write(f"Browsed {i} families\n")

        for fam2, homeo in fam_homeo.items():

            #ideally should check bijectivity but in practice rare that it's not, so this is faster
            if not fam1.isdisjoint(fam2):
                d_homeo[outgr] = homeo
                break
        i += 1

    return d_homeo


def load_combin(input_file, genes):

    """
    Loads family SCORPiOs family combination file (get family correspondance across multiple outgr)

    Args:
        input_file (str): path to the family combination file
        genes (set): list of the genes of the outgroup used as reference

    Returns:
        dict: genes in the reference outgroup to genes in the non-reference outgroups
    """
    combin = {}
    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):

            if i == 0:
                continue

            line = line.strip().split('\t')
            genes_all, best_graph = line[:-1], line[-1]

            if i == 1:
                for j, genes_ref in enumerate(genes_all):
                    if genes.intersection(set(genes_ref.split(','))):
                        break

            if genes.isdisjoint(set(best_graph.split(','))):
                genes_ref = genes_all[j].split(',')
                if genes_ref != []:
                    for gene in genes_ref:
                        combin[gene] = best_graph.split(',')
    return combin



def write_counter(counter, output_file):

    """
    Writes a dict object to file

    Args:
        counter (dict): input dict
        output_file (str): name of the output file
    """

    with open(output_file, 'w') as outfile:
        for key in counter:
            outfile.write(f'{counter[key]} {key}\n')


def write_output(d_homeo, ctrees, output_all, output_incons):

    """
    Write numbers of trees per homeologs

    Args:
        d_homeo (dict): family to homeologs correspondance
        ctrees (dict): family to synteny-sequence conflicts
        output_all (str): output file to write all considered families
        output_incons (str): output file to write sequence-synteny inconsistent families
    """

    tot = Counter(d_homeo.values())
    write_counter(tot, output_all)
    incons = {k:v for k, v in d_homeo.items() if ctrees[k] == "inconsistent"}
    incons = Counter(incons.values())
    write_counter(incons, output_incons)


if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--orthotable', nargs='+', help='SCORPiOs families file(s)',
                        required=True)

    PARSER.add_argument('-homeo', '--homeologs', help='Paralogy Map file (colored_ancGenes.tsv)',
                        required=True)

    PARSER.add_argument('-s', '--summary', help='SCORPiOs trees-synteny consistency summary',
                        required=True)

    PARSER.add_argument('-a', '--accepted', help='SCORPiOs accepted corrections',
                        required=True)

    PARSER.add_argument('-oa', '--out_all', required=False, default="out_all_trees.txt")

    PARSER.add_argument('-oi', '--out_incons', required=False, default="out_inconsistent_trees.txt")

    PARSER.add_argument('--is_outgroup', required=False, action='store_true')

    PARSER.add_argument('-f', '--fomt', help="genes coordinates file format",
                        required=False, default='bed')

    PARSER.add_argument('-c', '--combin', help="SCORPiOs family combination across outgroups.",
                        required=False, default=None)

    ARGS = vars(PARSER.parse_args())

    ACC = load_acc(ARGS["accepted"])
    CTREES = load_summary(ARGS["summary"], ACC)

    sys.stderr.write('Loading orthotable(s)...')
    D_OUTGR = {}
    for INPUT_FILE in ARGS["orthotable"]:
        D_OUTGR.update(load_outgr_fam(INPUT_FILE, CTREES))
    sys.stderr.write('ok\n')

    if not ARGS["is_outgroup"]:
        sys.stderr.write('Loading ancestral karyotype...')
        PM = load_pm(ARGS["homeologs"])
        sys.stderr.write('ok\n')

        sys.stderr.write("Transferring homeologs from the ancestral karyotype to SCORPiOs families "
                         f"({len(D_OUTGR)} families)...\n")

        HOMEOLOGS = outgroup_genes_to_homeologs(D_OUTGR, PM)
        sys.stderr.write('ok\n')

    else:
        GENOME = Genome(ARGS["homeologs"], ARGS["fomt"])
        COMBIN = None
        if ARGS["combin"] is not None:
            COMBIN = load_combin(ARGS["combin"], {g.names[0] for g in GENOME})
        HOMEOLOGS = {}
        GENOME = GENOME.genes_list
        for chrom in GENOME:
            for GENE in GENOME[chrom]:
                GENE = GENE.names[0]
                if GENE in CTREES:
                    HOMEOLOGS[GENE] = chrom
                elif COMBIN is not None and GENE in COMBIN:
                    for g in COMBIN[GENE]:
                        if g in CTREES:
                            HOMEOLOGS[g] = chrom

    write_output(HOMEOLOGS, CTREES, ARGS["out_all"], ARGS["out_incons"])
