#!/usr/bin/env python

"""
    Module with functions to load a genome from a .bed (or a .bz2 in DYOGEN format) gene file.
"""

import sys
import collections
import enum
import bz2


Gene = collections.namedtuple("Gene", ['chromosome', 'beginning', 'end', 'names'])
GenePosition = collections.namedtuple("GenePosition", ['chromosome', 'index'])


def is_bz2(filename):

    """
    Checks if file extension is bz2 (looks at file extension only, not its encoding,
    could be improved).

    Arg:
        filename (str): input file name

    Returns
        boolean: True if extension is bz2, False otherwise.
    """

    return filename.split('.')[-1] == 'bz2'


def toint(chr_name):

    """
    Converts the input to an integer, if possible. Otherwise leave the name unchanged,
    as str.

    Args:
        chr_name (str): String to convert, for instance a chromosome name.

    Returns:
        int or str: Converted input if possible, input otherwise
    """

    try:
        return int(chr_name)
    except TypeError:
        return None
    except ValueError:
        return sys.intern(chr_name) # for efficiency, same chr_names will be the same object


class ContigType(enum.Enum):

    """
    Enum grouping all possible values describing the type of a contig.
    """

    Chromosome = 'Chromosome'
    Mitochondrial = 'Mitochondrial'
    Scaffold = 'Scaffold'
    Random = 'Random'


def contig_type(chr_name):

    """
    Deduces the type of a contig from its name.

    Arg:
        chr_name (str): Name of the contig

    Returns:
        ContigType object: The type of the contig, either Chromosome, Mitochondrial, Scaffold or
                           Random
    """

    #if it can be converted to an integer it is usually a chromosome
    try:
        idx = int(chr_name)
        if idx < 100:
            return ContigType.Chromosome
        else:
            return ContigType.Scaffold

    except ValueError:

        chr_name_low = chr_name.lower()

        # chromosomes named that usually that corresponds to a scaffold that concatenates all
        # unassembled fragments
        if "rand" in chr_name_low or chr_name in ["UNKN", "Un", None]:
            return ContigType.Random

        # mitochondrial
        for mito_name in ['mt', 'mitochondrion']:
            if mito_name in chr_name_low:
                return ContigType.Mitochondrial

        #usual scaffold tags
        keys = ["cont", "scaff", "ultra", "reftig", "_", "un", "gl", "ki", "jh", "aaex", "aadn"]
        for tag in keys:
            if tag in chr_name_low:
                return ContigType.Scaffold
        if (chr_name in ["U", "E64", "2-micron"]) or chr_name.endswith("Het"):
            return ContigType.Scaffold

        else:
            # sex chromosomes
            return ContigType.Chromosome


class Genome:

    """
    Object representing genomic position of genes in a species, as loaded from a .bed (or 
    in DYOGEN format) gene file. Can load bzipped (.bz2) files.

    Attributes:
        name (str): name of the input gene file
        genes_list (dict): For each chromosome (key), a list of `Gene` namedtuples.
        chr_list (dict): For each ContigType (key), list of chromosomes with this type (value).
        dict_genes (dict): For each gene name (key), its position in a `GenePosition` namedtuple.
    """

    def __init__(self, fichier, file_format):

        """
        Inits a Genome Object from a given gene file.

        Arg:
            fichier (str): name of the genes coordinates file.
            format (str, optional): specify the format
        """

        sys.stderr.write("Loading genome of {} ...\n".format(fichier))

        if is_bz2(fichier):
            open_function = bz2.open
        else:
            open_function = open

        with open_function(fichier, 'rt') as infile:

            # list of genes per chromosome
            self.genes_list = collections.defaultdict(list)

            if file_format.upper() == "DYOGEN":

                # genesST or genes in DYOGEN format
                # CHR BEG END STRAND NAME
                # or CHR BEG END STRAND NAME ShortestTranscript

                for line in infile:
                    line = line.replace('\n', '').split('\t')
                    if len(line) == 5:
                        (chrom, beg, end, _, gene_names) = line
                    else:
                        assert len(line) == 6, "The genome file {} is not in the expected DYOGEN\
                                                .bz2 format, please check\n".format(fichier)

                        (chrom, beg, end, _, gene_names, _) = line

                    (beg, end) = (int(beg), int(end))
                    self.add_gene(gene_names.split(), chrom, beg, end)

            #if file is not bzip2, we assume it is a minimal .bed
            elif file_format.upper() == 'BED':

                # .bed minimal: "CHR BEG END NAMES"
                ###################################

                for line in infile:
                    line = line.replace('\n', '').split('\t')

                    assert len(line) == 4, "The genome file {} is not in the expected minimal .bed\
                                            format, please check\n".format(fichier)
                    (chrom, beg, end, gene_names) = line
                    self.add_gene(gene_names.split(), chrom, int(beg), int(end))

            else:
                raise ValueError("{} format is not supported\n".format(file_format))


        self.name = fichier

        self.init_other_attributes()


    def init_other_attributes(self):

        """
        Inits the genes and chromosomes dictionaries.
        """

        self.dict_genes = {}
        self.chr_list = collections.defaultdict(list)

        for chrom in self.genes_list:

            # associate a gene name to its location in the genome
            self.genes_list[chrom].sort() # sort by chromosome first and gene beg after
            for (idx, gene) in enumerate(self.genes_list[chrom]):
                for name in gene.names:
                    self.dict_genes[name] = GenePosition(chrom, idx)

            # classification into chromosomes/scaffolds/random contigs
            chr_t = contig_type(chrom)
            self.chr_list[chr_t].append(chrom)

        # sort chromosome names
        for chr_t in self.chr_list:
            self.chr_list[chr_t].sort(key=lambda v: (isinstance(v, str), v))


    def add_gene(self, names, chromosome, beg, end):

        """
        Adds a gene to the genes_list.

        Args:
            names (list): list of gene names
            chromosome (str): chromosome name
            beg, end (int): start and end positions of the gene
        """

        chromosome = toint(chromosome)
        self.genes_list[chromosome].append(Gene(chromosome, beg, end,
                                                tuple(sys.intern(name) for name in names)))


    # return all genes
    def __iter__(self):

        """
        Iterates over all genes, ordered by chromosome first and gene beg after.
        """

        for genes in self.genes_list.values():
            for gene in genes:
                yield gene


if __name__ == '__main__':
    sys.exit()
