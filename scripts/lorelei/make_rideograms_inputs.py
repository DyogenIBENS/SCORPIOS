#!/usr/bin/env python

"""
    Writes files that can be read by RIdeograms to draw karyotype with overlaid features.

    Example:

        $ python -m scripts.make_rideograms_inputs TODO
"""

import sys
import argparse

import roman

from scripts.synteny.mygenome import Genome


def strip_chr_name(chr_name):

    """
    Tries to strip a chr name so that it can be converted to int for RIdeograms.

    Args:
        chr_name (str): chr name

    Returns:
        (str): stripped chr name
    """

    if isinstance(chr_name, str) and "HiC" in chr_name:
        chr_name = chr_name.replace("HiC_scaffold_", "")

    elif isinstance(chr_name, str) and "chr" in chr_name:
        chr_name = chr_name.replace("chr", "")

    elif isinstance(chr_name, str) and "group" in chr_name:
        chr_name = chr_name.replace("group", "")

    elif isinstance(chr_name, str) and "LG" in chr_name:
        chr_name = chr_name.replace("LG", "")

    return chr_name

def make_karyo(genesfile, output, fomt='bed'):

    """
    Makes a karyotype file for drawing with RIdeograms from a bed file with genes coordinates.

    Args:
        genesfile (str): input file with genes coordinates
        output (str): output file name
        fomt (str, optional): input format .bed or dyogen format

    Returns:
        (scripts.synteny.mygenome.Genome): genome of the species for which to extract classes
        (list): ordered set of chromosomes.

    """

    genome = Genome(genesfile, fomt)
    dgenes = genome.genes_list
    karyo = []
    for chrom in dgenes:

        new_chrom = strip_chr_name(chrom)

        try:
            new_chrom = roman.fromRoman(new_chrom)
        except (TypeError, roman.InvalidRomanNumeralError):
            pass

        try:
            karyo.append((int(new_chrom), len(dgenes[chrom]), chrom))
        except ValueError:
            sys.stderr.write(f"Warning: chromosome {chrom} could not be converted to an integer,"
                             "it won't be drawn.\n")
            continue

    karyo = sorted(karyo)

    with open(output, 'w') as outfile:
        outfile.write("Chr\tStart\tEnd\n")

        for chrom in karyo:
            chrom, lg, _ = chrom
            outfile.write(f"{chrom}\t0\t{lg+1}\n")

    return genome, karyo


def load_features(genome, features_file, to_load=None):
    """
    Loads a 3-columns tab-delimited file with gene_family_name, genes and gene_family class.

    Args:
        genome (scripts.synteny.mygenome.Genome): genome of the species for which to extract classes
        features_file (str): path to the input file with gene family classes
        to_load (list, optional): load only genes of given classes

    Returns:
        dict: for genes in the input genome (key) gives the gene family class (value)
    """
    feat = {}
    sp_genes = {g.names[0] for g in genome}
    with open(features_file, 'r') as infile:

        for line in infile:

            _, descendants, classif = line.strip().split('\t')


            if to_load is None or classif in to_load:

                try:
                    classif = int(classif)
                except ValueError:
                    if classif == "Inconsistent":
                        classif = 1
                    else:
                        pass

                genes = set(descendants.split()).intersection(sp_genes)

                if genes:
                    genes = list(genes)
                    for gene in genes:
                        feat[gene] = classif
    return feat


def features_to_ide(genome, features_file, karyo, output, to_load=None):

    """
    Writes the gene to gene family class to file, to use as input to RIdeogram.

    Args:
        genome (scripts.synteny.mygenome.Genome): genome of the species for which to extract classes
        features_file (str): path to the input file with gene family classes
        karyo (list): ordered set of chromosome
        output (str): name for the output file
        to_load (list, optional): load only genes of given classes

    """

    feat = load_features(genome, features_file, to_load)

    dgenes = genome.genes_list

    with open(output, 'w') as out:
        out.write("Chr\tStart\tEnd\tValue\n")

        for k in karyo:
            chrom_int, _, chrom = k

            for i, gene in enumerate(dgenes[chrom]):
                name = gene.names[0]

                if name in feat:
                    classif = feat[name]
                    start = i
                    stop = i + 1

                    out.write(f"{chrom_int}\t{start}\t{stop}\t{classif}\n")


if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--input_incons', help='', required=True)

    PARSER.add_argument('-g', '--genesfile', help='', required=True)

    PARSER.add_argument('-f', '--format', help='', required=False, default="bed")

    PARSER.add_argument('-k', '--outfile_karyo', help='Output file rideogram karyotype.',
                        required=False, default="out")

    PARSER.add_argument('-o', '--outfile_features', help='Output file rideogram overlaid data.',
                        required=False, default="out")

    PARSER.add_argument('-t', '--to_load', nargs='*', required=False, default=None)

    ARGS = vars(PARSER.parse_args())

    GENOME, KARYO = make_karyo(ARGS["genesfile"], ARGS["outfile_karyo"], fomt=ARGS["format"])

    features_to_ide(GENOME, ARGS["input_incons"], KARYO, ARGS["outfile_features"], ARGS["to_load"])
