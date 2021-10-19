#!/usr/bin/env python

"""
    Module with functions to compare duplicated segments in 2 duplicated species.
"""

import sys
from collections import defaultdict
import itertools

import numpy as np

from . import utilities as ut


def to_dup_segments(fams):

    """
    Transforms a list of GeneFamilies (i.e entries of a given duplicated species in a window of the
    OrthologyTable) into a DupSegments object.

    A DupSegments object consist in:

    - (i): a list of names of each gene family, given by the corresponding outgroup gene in the
      OrthologyTable
    - (ii): a list of all genomic segments with a gene copy in the duplicated species
    - (iii): a binary matrix, representing absence/presence of a duplicated gene copy in
      each genomic segment. Columns are genomic segments, with order given in list (ii). Rows are
      gene families.
    - (iv): a dictionary, giving for each '1' in the matrix, corresponding duplicated species
      gene names
    - (v): a list keeping track of discarded families segments threadings

    Arg:
        fams (list of GeneFamily objects): object to transform

    Returns:
        DupSegments: the transformed object
    """

    #(i) extract gene name in the outgroup family ids
    family_ids = ut.get_all_outgr_names(fams)

    #(ii) split chromosome that have ohnologs in two different regions
    _, splits = ut.split_chr_with_ohnologs(fams)

    #(ii) get a list of all genomic segments with a gene in the duplicated species
    chromosomes = [i for i in ut.get_all_chromosomes_involved(fams) if i not in splits.keys()]

    #(ii) add split chromosomes to the list of genomic segments
    if splits:
        chromosomes += [''.join(map(str, x)) for x in itertools.product(sorted(splits.keys()),\
                                                                        ['_a', '_b'])]

    #(iii) build a binary matrix representing absence/presence of a duplicated gene in each segment
    #First init the matrix
    n_rows = len(fams)
    n_col = len(chromosomes)
    matrix = np.zeros((n_rows, n_col), dtype=np.bool)
    genes = defaultdict(dict)

    #Then fill the matrix (iii) and the corresponding dict (iv)
    for i in range(n_rows):
        for gene in fams[i].all_duplicate_genes:
            chrom = gene.chromosome
            if chrom:

                if chrom in splits:
                    chrom += '_'+splits[chrom][gene.name]

                j = chromosomes.index(chrom)
                matrix[i, j] = 1
                genes[i][j] = genes[i].get(j, []) + [gene.name]

    #init the DupSegments object
    dup_seg = DupSegments(family_ids, chromosomes, matrix, genes)

    #sort to have segments with most genes in the first columns
    dup_seg.sort()

    return dup_seg


class DupSegments:

    """
    Object to represent a list of GeneFamilies (i.e entries of a given duplicated species in a
    window of the OrthologyTable). This object is used to perform duplicated segment threading and
    synteny similarity comparisons in pairs of duplicated species.

    Attributes:
        family_ids (list of str): Names of gene families, given by the outgroup gene in the
                                  OrthologyTable

        chromosomes (list of str): Names of genomic segments with a gene copy in the duplicated
                                   species

        matrix (numpy.array): A binary matrix, representing absence/presence of a duplicated gene
                              copy in each genomic segment. Columns are genomic segments, with
                              order given in list (ii). Rows are gene families.

        genes_dict (dict): For each '1' in the matrix, the corresponding duplicated species gene
                           name(s)
    """

    def __init__(self, family_ids, chromosomes, matrix, genes_dict):

        """
        Inits DupSegments with its pre-computed attributes value.
        """

        self.family_ids = family_ids #id of each gene family (name of the outgroup gene)
        self.chromosomes = chromosomes #all genomic segments of duplicated genes
        self.matrix = matrix #matrix describing ortholog genes presence/absence
        self.genes_dict = genes_dict #corresponding dictionary with gene names

        self.discard = [] #to store where a duplicated gene was discarded during threading


    def sort(self):
        """
        Orders duplicated segments by descending number of genes.
        """

        index_map = {}
        nb_genes = np.count_nonzero(self.matrix, axis=0)

        #sort using mergesort to have a deterministic sort
        sorted_ind = np.argsort(nb_genes, kind='mergesort')[::-1]
        for new_ind, prev_ind in enumerate(sorted_ind):
            index_map[prev_ind] = new_ind

        new_genes_d = defaultdict(dict)

        for keyx in self.genes_dict.keys():

            for keyy in self.genes_dict[keyx].keys():

                new_genes_d[keyx][index_map[keyy]] = self.genes_dict[keyx][keyy]

        self.chromosomes = [self.chromosomes[i] for i in sorted_ind]
        self.matrix = self.matrix[:, sorted_ind]
        self.genes_dict = new_genes_d


    def all_reduce_in_two_blocks(self):
        """
        Gives a list of all possible ways to thread duplicated segments together.

        Returns:
            nested list: all possible threadings.

            For instance, a threading given by [[0, 2], [1, 3]] means segments 0 and 2 threaded
            together to form an ancestrally duplicated region track1 and segment 1 and 3 threaded
            together to form track2.
        """

        #total number of segments
        n_col = len(self.matrix[0, :])
        j = 0

        if n_col == 1:
            return [[[0], []]]

        if n_col == 2:
            return [[[0], [1]]]

        #start by placing segment with max number of genes in each of the 2 tracks
        already_placed = 2
        threadings = [[[0], [1]]]

        #continue until we have no more unthreaded segment
        while n_col - j >= already_placed + 1:

            current_threadings = []

            #browse all already computed segment threading possibilities
            for threading in threadings:

                #segment threaded together in each of the two tracks
                seg_track1 = threading[0]
                seg_track2 = threading[1]

                #get positions with a gene copy in each track and in current segment
                major1, major2, curr_seg = np.count_nonzero(self.matrix[:, seg_track1], axis=1),\
                                           np.count_nonzero(self.matrix[:, seg_track2], axis=1),\
                                           self.matrix[:, j+already_placed]

                #threading with segments in track1 possible if no ohnolog gene clash
                if np.all(major1+curr_seg < 2):

                    current_threadings += [[seg_track1+[j+already_placed], seg_track2]]

                #threading with segments in track2 possible if no ohnolog gene clash
                if np.all(major2+curr_seg < 2):

                    current_threadings += [[seg_track1, seg_track2+[j+already_placed]]]


            #update all threadings if we successfully threaded current segment j+2
            if current_threadings != []:

                threadings = current_threadings

            j += 1

        return threadings


    def update_discard(self, threading):

        """
        Updates the discard attribute.

        Args:
            threading (nested list): threading scenario.
        """

        #for all gene families in the window
        for loc in self.genes_dict.keys():

            #for all associated segments
            for chrom in self.genes_dict[loc].keys():

                #if segment could not be threaded
                if chrom not in threading[0] and chrom not in threading[1]:

                    #store gene family position in the window
                    if self.genes_dict[loc][chrom]:
                        self.discard.append(loc)


    def update_orthologies(self, dup_seg_sp2, score, threading2sp, all_orthologies):

        """
        Stores genes in identified orthologous duplicated segment. Fills `all_orthologies`
        in-place.

        Args:
            dup_seg_sp2 (DupSegments): corresponding duplicated segments in species 2
            score (float): delta score of synteny similarity (diff. of the 2 orthology scenarios)
            threading2sp (nested list): duplicated segment threading for each species
            all_orthologies (dict): stores orthologies, for each family (key) gives a tuple (value)
                                    with the confidence score and predicted orthologs.
        """

        #group of orthologous segments based on score sign (deltaS)
        pairs = [(0, 0), (1, 1)]
        if score > 0:
            pairs = [(0, 1), (1, 0)]

        #update discard attribute
        threads_sp1, threads_sp2 = threading2sp
        self.update_discard(threads_sp1)
        dup_seg_sp2.update_discard(threads_sp2)

        for loc in self.genes_dict:

            family_id = self.family_ids[loc]


            #update orthology for family if deltaS score is higher
            if family_id not in all_orthologies or abs(score) >= all_orthologies[family_id][0]:

                #ignore family where a gene was discarded
                if loc not in self.discard and loc not in dup_seg_sp2.discard:

                    new_ortho = []
                    sp1_nogene = 0
                    sp2_nogene = 0

                    #for matched orthologous segments
                    for pair in pairs:

                        genes1 = []
                        genes2 = []

                        #gene(s) in sp1 in orthologous segment
                        for chrom in threads_sp1[pair[0]]:

                            genes1 += self.genes_dict[loc].get(chrom, [])

                        if not genes1:

                            sp1_nogene += 1

                        #gene(s) in sp2 in orthologous segment
                        for chrom in threads_sp2[pair[1]]:

                            genes2 += dup_seg_sp2.genes_dict[loc].get(chrom, [])

                        if not genes2:

                            sp2_nogene += 1

                        new_ortho.append((genes1, genes2))

                    #check that at least 1 gene in each species for this family and update
                    if sp2_nogene < 2 and sp1_nogene < 2:

                        all_orthologies[family_id] = [abs(score), new_ortho]


    def get_score(self, dup_seg_sp2, tree_orthos, threadingsp1, threadingsp2):

        """
        Computes the two delta scores between tracks of threaded duplicated segments in 2 species.

        Args:
            dup_seg_sp2 (DupSegments): Corresponding duplicated segments in species 2
            tree_orthos (dict): Orthologous gene pairs in sp1 and sp2, defined from molecular
                                evolution
            threadingsp1 (nested list): duplicated segment threading for species 1
            threadingsp2 (nested list): duplicated segment threading for each species 2
        Returns:
            tuple: tuple of 2 floats, delta score based on the 'pattern of retentions and losses'
            and delta score based on 'syntenic neighbours'

        """

        pattern_score = self.retention_loss_score(dup_seg_sp2, threadingsp1, threadingsp2)

        ortho_score = self.orthologs_score(dup_seg_sp2, tree_orthos, threadingsp1, threadingsp2)

        return (pattern_score, ortho_score)


    def retention_loss_score(self, dup_seg_sp2, threadingsp1, threadingsp2):

        """
        Computes delta score based on the 'pattern of retentions and losses' between tracks of
        threaded duplicated segments in 2 species.


        Args:
            dup_seg_sp2 (DupSegments): Corresponding duplicated segments in species 2
            threadingsp1, threadingsp2 (nested list): duplicated segments threading
                                                      for each species

        Returns:
            float: delta score based on the 'pattern of retentions and losses'
        """

        ##sp1##

        #binary vector representing absence/presence of genes in track1 threaded segments
        e1_a = np.bitwise_or.reduce(self.matrix[:, threadingsp1[0]], axis=1)
        length = len(e1_a)

        #binary vector representing absence/presence of genes in track2 threaded segments
        if self.matrix.shape[1] > 1:
            e1_b = np.bitwise_or.reduce(self.matrix[:, threadingsp1[1]], axis=1)
        else:
            e1_b = np.zeros(length, dtype=np.bool)

        ##sp2##

        #binary vector representing absence/presence of genes in track1 threaded segments
        e2_a = np.bitwise_or.reduce(dup_seg_sp2.matrix[:, threadingsp2[0]], axis=1)

        #binary vector representing absence/presence of genes in track2 threaded segments
        if dup_seg_sp2.matrix.shape[1] > 1:
            e2_b = np.bitwise_or.reduce(dup_seg_sp2.matrix[:, threadingsp2[1]], axis=1)
        else:
            e2_b = np.zeros(length, dtype=np.bool)

        #Number of common retentions and losses in scenario1
        #computed with 1 - bit-wise XOR for efficiency
        matches_aa = len(np.bitwise_xor(e1_a, e2_a).nonzero()[0])
        matches_bb = len(np.bitwise_xor(e1_b, e2_b).nonzero()[0])
        sim_scenario1 = 1 - ((matches_aa + matches_bb) / float(length*2))

        #Number of common retentions and losses in scenario2
        matches_ab = len(np.bitwise_xor(e1_a, e2_b).nonzero()[0])
        matches_ba = len(np.bitwise_xor(e1_b, e2_a).nonzero()[0])
        sim_scenario2 = 1 - ((matches_ab + matches_ba) / float(length*2))

        #delta
        delta_score = sim_scenario2 - sim_scenario1

        return delta_score


    def orthologs_score(self, dup_seg_sp2, tree_orthos, threadingsp1, threadingsp2):

        """
        Computes delta score based on 'syntenic neighbours' between tracks of threaded duplicated
        segments in 2 species.


        Args:
            dup_seg_sp2 (DupSegments): Corresponding duplicated segments in species 2
            tree_orthos (dict): Orthologous gene pairs in sp1 and sp2, defined from molecular
                                evolution
            threadingsp1, threadingsp2 (nested list): duplicated segments threading for each
                                                      species

        Returns:
            float: delta score based on 'syntenic neighbours'
        """

        matches_aa_bb, matches_ab_ba = 0, 0

        #browse genefamilies with a gene copy in species 1
        for loc in self.genes_dict:

            #browse corresponding genomic segments
            for chrom in self.genes_dict[loc]:

                genes1 = self.genes_dict[loc][chrom]

                donea = False
                doneb = False

                #browse gene copies on this genomic segment for this gene family
                for gene1 in genes1:

                    #continue if it is in pre-computed orthologies
                    if gene1 in tree_orthos:

                        #its pre-computed ortholog(s) in species 2
                        pred = tree_orthos[gene1]

                        #matched orthologous segments in species 2
                        chroms_sens1, chroms_sens2 = threadingsp2
                        if chrom in threadingsp1[1]:
                            chroms_sens2, chroms_sens1 = threadingsp2

                        #see if genes matched in scenario 1 are orthologous
                        if not donea:

                            found = check_orthology(chroms_sens1, dup_seg_sp2, pred, loc)

                            if found:
                                matches_aa_bb += 1
                                donea = True

                        #see if genes matched in scenario 2 are orthologous
                        if not doneb:

                            found = check_orthology(chroms_sens2, dup_seg_sp2, pred, loc)

                            if found:
                                matches_ab_ba += 1
                                doneb = True


        #normalise counts and compute delta
        tot = len(self.matrix[:, 0])*2

        scenario1 = (matches_aa_bb) / float(tot)
        scenario2 = (matches_ab_ba) / float(tot)

        delta_score = scenario2 - scenario1

        return delta_score


def check_orthology(orthologous_chroms, dup_seg_sp2, ortho_genes, loc):

    """
    Checks if there is a pre-computed orthology relation between genes of matched duplicated
    segments in 2 species for the family `loc`.

    Args:
        orthologous_chroms (list): list of orthologous segments in species 2
        dup_seg_sp2 (DupSegments): duplicated segments object for species 2
        ortho_genes (list): list of orthologous genes in species 2 for a gene in species 1
        loc (int): index of the gene family

    Returns:
        bool: True if there is a pre-computed gene orthology, False otherwise.
    """

    found = False
    i = 0

    #browse all pre-computed orthologs
    while i < len(ortho_genes) and not found:

        #current pre-computed ortohlog
        ortho = ortho_genes[i]

        #search if current ortholog is on matched orthologous segments
        for chrom_sp2 in orthologous_chroms:

            if loc in dup_seg_sp2.genes_dict and chrom_sp2 in dup_seg_sp2.genes_dict[loc]:

                if ortho in dup_seg_sp2.genes_dict[loc][chrom_sp2]:

                    found = True

        i += 1

    return found

if __name__ == '__main__':
    sys.exit()
