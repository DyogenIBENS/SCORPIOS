#!/usr/bin/env python

"""
    Module with functions to load and write a duplicated ingroups-outgroup orthology table.
"""

import sys

import collections
import numpy as np


GeneSpeciesPosition = collections.namedtuple("GenePosition", ['name', 'chromosome', 'index'])


def outgr_chromosomes(chr_file):

    """
    Reads a simple file with a single entry on each line (for instance chrom names) on each line.

    Args:
        chr_file (str): input file name

    Returns:
        list: entry on each line of the file
    """

    chroms = []
    with open(chr_file, 'r') as infile:
        for line in infile:
            chroms.append(line[:-1])
    return chroms


def complete_load_orthotable(table_file, chrom_outgr, species, load_no_position_genes=False):

    """
    Loads entries for one duplicated species `species` in the orthologytable, corresponding to
    chromosome `chrom_outgr`in the outgroup, as a list of `GeneFamily` objects.

    Args:
        table_file (str): Name of the orthologytable file
        chrom_outgr (str): Name of the considered outgroup chromosome
        species (str): Name of the considered duplicated species

    Returns:
        list of `GeneFamily` objects
    """

    tmp = []
    list_of_genefam = []

    with open(table_file, 'r') as infile:

        for i, line in enumerate(infile):

            #Read header to get species in each column
            if i == 0:
                first_line = line
                spl = first_line.strip().split('\t')
                ind = spl.index(species)
                i = 1

            #Add other lines
            elif line.split('\t')[0] == chrom_outgr:

                spl = line[:-1].split('\t')
                outgr_chr, outgr_position, outgr_gname = spl[:3]
                all_genes = []
                tmp_chrom = []

                #Extract infos for all genes in duplicated species `species`
                for gene in spl[ind].split('/'):

                    if gene:
                        name, info = gene.split('|')

                        if info != '-':

                            chrom, position = info.split('-')
                            all_genes.append(GeneSpeciesPosition(name, chrom, int(position)))

                            if chrom not in tmp_chrom:
                                tmp_chrom.append(chrom)

                        elif load_no_position_genes:
                            all_genes.append(GeneSpeciesPosition(name, '', ''))


                #Add gene family
                tmp = GeneFamily(outgr_gname, outgr_chr, outgr_position, all_genes, tmp_chrom)
                list_of_genefam.append(tmp)

    return list_of_genefam


def load_orthotable(table_file):

    """
    Simplified loading function for the orthologytable, in order to get outgroup genes in the
    orthology table and all duplicated species gene copies in the corresponding family.

    Args:
        table_file (str): Input file name.

    Returns:
        orthotable (dict): Correspondence between genes in the outgroup (key) and duplicated
        species genes in its family (value).
    """

    orthotable = {}

    with open(table_file, 'r') as infile:

        for i, line in enumerate(infile):

            #skip header
            if i > 0:
                spl = line[:-1].split('\t')
                outgr_gene = spl[2]
                orthotable[outgr_gene] = []

                for species in spl[3:]:
                    genes = species.split('/')

                    if genes != ['']:
                        for gene in genes:
                            name, _ = gene.split('|')
                            orthotable[outgr_gene].append(name)
    return orthotable


def light_load_orthotable(table_file):

    """
    Another simplified loading function for the orthologytable, in order to only get outgroup
    genes in the orthology table and their corresponding chromosomes.

    Args:
        input_file (str): Input file name.

    Returns:
        names (dict): Correspondence between chromosome of the outrgoup (key) and its genes in the
        orthology table (value). Genes are given in order of along each chromosome.
    """

    names = {}

    with open(table_file, 'r') as infile:

        for i, line in enumerate(infile):

            #skip header
            if i != 0:

                spl = line[:-1].split('\t')
                outgr_gname = spl[2]
                outgr_chr = spl[0]
                names[outgr_chr] = names.get(outgr_chr, [])
                names[outgr_chr].append(outgr_gname)

    return names


class GeneFamily:

    """
    Stores an entry in the orthology table for one duplicated species.

    Attributes:
        outgr_genename (str): name of the outgroup gene, giving an unique IDs to the family
        outgr_chr (str): name of the chromosome of the outgroup gene
        outgr_position (int): index of the outgroup gene on its chromosome
        all_duplicate_genes (list of `GeneSpeciesPosition`): gene copies in the duplicated
                            species and their genomic location
        involved_chromosomes (list of str): list of chromosomes in the duplicated species with a
                             gene copy

    Note:
        No public method, used as a structure to store data. `GeneFamily` objects are manipulated
        in lists with functions on `GeneFamily`lists defined below for better readability of
        manipulations.

    """

    def __init__(self, outgr_genename, outgr_chr, outgr_position, all_duplicate_genes,
                 involved_chromosomes):
        """
        Inits a `GeneFamily` object
        """
        self.outgr_genename = outgr_genename
        self.outgr_chr = outgr_chr
        self.outgr_position = int(outgr_position)
        self.all_duplicate_genes = all_duplicate_genes
        self.involved_chromosomes = involved_chromosomes


def get_all_chromosomes_involved(list_of_genefam):

    """
    Gets all chromosome with a gene copy in a list of `GeneFamily` objects.

    Args:
        list_of_genefam (list of `GeneFamily`): input list of `GeneFamily`

    Returns:
        list of str: list of chromosome names
    """

    chrom = [j.chromosome for i in list_of_genefam for j in i.all_duplicate_genes\
             if j.chromosome]
    used = set()
    return [x for x in chrom if x not in used and (used.add(x) or True)]


def get_all_outgr_pos(list_of_genefam):

    """
    Gets chromosomal location index of all outgroup genes in a list of `GeneFamily` objects.

    Args:
        list_of_genefam (list of `GeneFamily`): input list of `GeneFamily`

    Returns:
        list of int: list of chromosomal indexes
    """

    return [i.outgr_position for i in list_of_genefam]


def get_all_outgr_names(list_of_genefam):

    """
    Gets gene names of all outgroup genes in a list of `GeneFamily` objects.

    Args:
        list_of_genefam (list of `GeneFamily`): input list of `GeneFamily`

    Returns:
        list of str: list of gene names
    """
    return [i.outgr_genename for i in list_of_genefam]


def get_all_chromosome_and_position(list_of_genefam):

    """
    Gets chromosome and chromosomal location index of all the duplicated species genes in a list
    of `GeneFamily` objects.

    Args:
        list_of_genefam (list of `GeneFamily`): input list of `GeneFamily`

    Returns:
        dict: for each chromosome (key), list of gene positions (value)
    """

    locations = {}
    loc = [(j.chromosome, j.index) for i in list_of_genefam for j in i.all_duplicate_genes\
           if j.chromosome]
    for chrom, pos in loc:
        locations.setdefault(chrom, []).append(pos)
    return locations


def split_chr_with_ohnologs(list_of_genefam):

    """
    Splits the duplicated species chromosomes in two separate regions if there are two ohnolgs on
    the same chromosome but more than 100 genes apart. This can potentially be the result of
    different duplicated chromosomes that fused together.

    Args:
        list_of_genefam (list of `GeneFamily`): input list of `GeneFamily`

    Returns:
        store (store): list tuples with an historic of ohnologs 100 genes apart on the same
                       chromosome
        splits (dict of dict): for each split chromosome (key1), each gene copy on it (key2) and
                               its corresponding after-split region ('a' or 'b')
    """

    clash = {}

    store = []
    splits = {}

    #Find self-clash (2 ohnologs) chr and extract most-distant position
    for i in list_of_genefam:
        seen = {}
        for j in i.all_duplicate_genes:
            if j.chromosome in seen:
                diff = np.array([abs(i[0] - j.index) for i in seen[j.chromosome]])
                maxi, argmax = np.max(diff), np.argmax(diff)
                if (j.chromosome not in clash and maxi > 100) or (j.chromosome in clash\
                                                             and maxi > clash[j.chromosome][0]):

                    clash[j.chromosome] = (maxi, (seen[j.chromosome][argmax][0], j.index))
                    store.append((maxi, seen[j.chromosome][argmax][1], j.name))

            else:
                seen[j.chromosome] = []

            seen[j.chromosome].append((j.index, j.name))

    if clash:

        #Classify in 2 groups based on their location on chromosome
        for i in list_of_genefam:
            for j in i.all_duplicate_genes:
                if j.chromosome in clash:
                    dist_a = abs(j.index - clash[j.chromosome][1][0])
                    dist_b = abs(j.index - clash[j.chromosome][1][1])
                    if dist_a < dist_b:
                        if j.chromosome not in splits:
                            splits[j.chromosome] = {}

                        splits[j.chromosome][j.name] = 'a'
                    else:
                        if j.chromosome not in splits:
                            splits[j.chromosome] = {}
                        splits[j.chromosome][j.name] = 'b'

    return store, splits


def insert_outgr_gene(list_of_genefam, ind, gene):

    """
    Inserts an outgroup gene in the orthology table (i.e in the list of `GeneFamily`).

    Args:
        list_of_genefam (list of `GeneFamily`): input list of `GeneFamily`
        ind (int): index to insert the gene in the list
        gene (`GeneSpeciesPosition` namedtuple): gene to add
    """

    to_add = GeneFamily(gene.name, gene.chromosome, gene.index, [], [])

    if ind == 0:
        new_table = [to_add]+list_of_genefam

    elif ind > len(list_of_genefam):
        new_table = list_of_genefam + [to_add]

    else:
        new_table = list_of_genefam[0:ind]+[to_add]+list_of_genefam[ind:]

    return new_table


def add_gene(list_of_genefam, ind, gene):

    """
    Adds a gene copy member of the duplicated species in the orthology table (i.e in the
    corresponding `GeneFamily`).

    Args:
        list_of_genefam (list of `GeneFamily`): input list of `GeneFamily`
        ind (int): family index to add the gene in the list
        gene (`GeneSpeciesPosition` namedtuple): gene to add
    """

    list_of_genefam[ind].all_duplicate_genes.append(gene)

    if gene.chromosome not in list_of_genefam[ind].involved_chromosomes:
        list_of_genefam[ind].involved_chromosomes.append(gene.chromosome)

    return list_of_genefam


def update_orthologytable(all_genefam, res_dict, sp_list):

    """
    Updates the orthology table by adding all newly found orthologies between ingroups and the
    outgroup. Inserts the new family in the orthology table (i.e in the list of `GeneFamily`).

    Args:
        all_genefam (nested dict): Stores the full orthology table. For each chromosome of the
                                   outgroup (key1), for each duplicated species (key2), a list of
                                   `GeneFamily` objects (value).
        res_dict (nested dict): Stores new orthologies. For each gene family (key1; represented by
                                its family id, the outgroup gene name), for each duplicated
                                species (key2), the corresponding gene copies in the duplicated
                                species (value).
        sp_list (list of str): list of duplicated species
    """

    #for each new family
    for outgr_gene in res_dict:

        chrom = outgr_gene.chromosome
        idx = outgr_gene.index

        if chrom in all_genefam: #should always be the case

            #for each duplicated species with a gene copy in the family
            for spec in res_dict[outgr_gene]:

                for gene in res_dict[outgr_gene][spec]:

                    #find index were to insert the family in the list
                    pos = get_all_outgr_pos(all_genefam[chrom][spec])
                    ind = find_closest(idx, pos, index=True)

                    #if not already an entry for this gene family
                    if idx not in pos:

                        #Insert outgr gene
                        all_genefam[chrom][spec] = insert_outgr_gene(all_genefam[chrom][spec],
                                                                     ind, outgr_gene)

                    #add current gene
                    all_genefam[chrom][spec] = add_gene(all_genefam[chrom][spec], ind, gene)


            #Also insert the `GeneFamily` for species without orthologs
            for spec in sp_list:

                if spec not in res_dict[outgr_gene]:

                    pos = get_all_outgr_pos(all_genefam[chrom][spec])

                    if idx not in pos:
                        all_genefam[chrom][spec] = insert_outgr_gene(all_genefam[chrom][spec],
                                                                     ind, outgr_gene)




def write_updated_orthotable(all_genefam, outgr, sp_list, chr_outgr, out, wsize=0,
                             filt_genes=None):

    """
    Writes an orthology table file from data stored in `all_genefam`.

    Args:
        all_genefam (nested dict): Stores the full orthology table. For each chromosome of the
                                   outgroup (key1), for each duplicated species (key2), a list of
                                   `GeneFamily` objects (value).
        outgr (str): name of the outgroup species
        sp_list (list of str): list of duplicated species
        chr_outgr (str): chromosome of the outgroup
        out (str): name of the output file
    """

    stats = {}

    with open(out, 'w') as outfile:

        #header
        outfile.write(outgr+'\t\t\t\t'+('\t').join(sp_list)+'\n')

        sp1 = sp_list[0]

        #for each outgroup chromosome
        for chrom in chr_outgr:
            length = len(all_genefam[chrom][sp1])

            #for each gene family
            for i in range(length):

                stats['families'] = stats.get('families', 0) + 1

                sp_with_genes = 0

                #outgroup gene information
                outgr_info = all_genefam[chrom][sp1][i]
                outgr_chr = outgr_info.outgr_chr
                outgr_genename = outgr_info.outgr_genename

                outgr_pos = outgr_info.outgr_position
                tmp = []
                tot_genes = 0
                #orthologs in all species
                for species in sp_list:

                    stats[species] = stats.get(species, 0)

                    genes = all_genefam[chrom][species][i].all_duplicate_genes

                    if not genes:
                        tmp.append('')

                    else:
                        content = []
                        for gene in genes:
                            stats[species] += 1
                            tot_genes += 1
                            if filt_genes and gene.name in filt_genes:
                                content.append(gene.name+'|-')
                            else:
                                content.append(gene.name+'|'+gene.chromosome+'-'+str(gene.index))
                        str_content = ('/').join(content)
                        tmp.append(str_content)

                        sp_with_genes += 1

                    teleost_genes = ('\t').join(tmp)

                outfile.write(outgr_chr+'\t'+str(outgr_pos)+'\t'+outgr_genename+'\t'+'\t')
                outfile.write(teleost_genes+'\n')

                if sp_with_genes == 1 or tot_genes == 2:
                    stats['families_1g'] = stats.get('families_1g', [])
                    stats['families_1g'].append(outgr_genename)

                elif wsize and length < wsize:
                    stats['small_chr'] = stats.get('small_chr', [])
                    stats['small_chr'].append(outgr_genename)

    return stats


def find_closest(number, number_list, index=False):

    """
    Finds, in a list of int `number_list`, the closest integer to `number`, or its index.
    Assumes the list is sorted. If two values are equally close to number, gives the smallest.

    Args:
        number (int): the input number to search
        number_list (list): the list of int to mine
        index (bool, optional): Whether the index of the closest element should be returned instead
                             of its value.

    Returns:
        int: closest number in list (or ist index if index is True)
    """

    closest = min(number_list, key=lambda x: abs(x-number))

    if index:
        ind = number_list.index(closest)
        if closest == number:
            return ind
        if closest < number:
            return ind+1
        return ind

    else:
        return closest

if __name__ == '__main__':
    sys.exit()
