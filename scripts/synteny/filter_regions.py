#!/usr/bin/env python

"""
    Module with functions to extract gene families having updated synteny information in SCORPiOs
    iteration n versus iteration n-1.
"""

import os
import sys

from .utilities import light_load_orthotable, find_closest


def get_modified_families(orthotable, orthotable_prev, corrected_fam, mapping_fam=None):

    """
    For OrthologyTables of two successive SCORPiOs iterations, find families with updated
    homologies in iteration n compared to iteration n-1.

    Updated families are either (i) a corrected tree or (ii) an outgroup gene in iteration n
    without duplicated species orthologs in iteration n-1.

    Args:
        orthotable (dict): gene families at iteration n
        orthotable_prev (dict): gene families at iteration n-1
        corrected_fam (list): list of corrected families
        mapping_fam (dict, optional): when multiple outgroups are used, a dictionary with
                                      correspondence of families ids across outgroups, useful when
                                      a tree was corrected using an other outgroup

    Returns:
        dict: modified gene families
    """

    modified_fam = {}
    for chrom in orthotable:
        modified_fam[chrom] = []
        for gene in orthotable[chrom]:
            if gene not in orthotable_prev[chrom]:
                modified_fam[chrom].append(gene)
            elif gene in corrected_fam:
                modified_fam[chrom].append(gene)
            elif mapping_fam and gene in mapping_fam:
                if set(mapping_fam[gene]).intersection(set(corrected_fam)):
                    modified_fam[chrom].append(gene)
    return modified_fam


def get_genes_to_keep(orthotable, modified_fam, windowsize):

    """
    Extracts all families with updated synteny information after a SCORPiOs iteration, i.e. all
    families within the same windows as a modified family.

    Args:
        orthotable (dict): gene families at iteration n
        modified_fam (dict): modified gene families

    Returns:
        dict: for each chromosome, families with updated synteny information
        list: flat list of updated families (outgroup gene name)
    """

    to_save = {}
    all_kept = []
    for chrom in orthotable:
        to_save[chrom] = []
        idx = [orthotable[chrom].index(i) for i in modified_fam[chrom]]
        for i, gene in enumerate(orthotable[chrom]):

            found = False
            if gene in modified_fam[chrom]:
                to_save[chrom].append((gene, i))
                found = True
                all_kept.append(gene)

            elif idx:
                closest = find_closest(i, idx)
                if abs(closest - i) <= windowsize:
                    to_save[chrom].append((gene, i))
                    found = True
                    all_kept.append(gene)

            if not found:
                to_save[chrom].append('-')

    return to_save, all_kept


def write_regions_file(fam_to_keep, outfile):

    """
    Writes a file with families that have updated synteny information after a SCORPiOs iteration.

    Args:
        fam_to_keep (dict): for each outgroup chromosome, families with updated synteny info.
        outfile (str): name of the output file.
    """

    with open(outfile, 'w') as out:
        for chrom in fam_to_keep:
            prev = ''

            for gene in fam_to_keep[chrom]:

                if (prev == '-' or not prev) and gene != '-':
                    out.write('START\n')

                if len(gene) > 1:
                    name, i = gene
                    out.write(chrom+'\t'+name+'\t'+str(i)+'\n')

                prev = gene


def read_combin_file(file_combin_graphs):

    """
    Reads the summary of graphs combination across outgroups. Corrected subtrees with another
    outgroup should be marked as an updated family for all outgroups.

    Args:
        file_combin_graphs (str): input summary file

    Returns:
        dict: for each gene in the current outgroup, the corresponding selected graph if from
              another outgroup
    """

    mapping = {}
    with open(file_combin_graphs, 'r') as infile:
        for i, line in enumerate(infile):
            if i != 0:
                line = line.strip().split('\t')
                selected_graphs = line[-1].split(',')
                for spec in line[:-1]:
                    for gene in spec.split(','):
                        if gene not in selected_graphs:
                            for selected_graph in selected_graphs:
                                mapping[gene] = mapping.get(gene, [])
                                mapping[gene].append(selected_graph)
    return mapping


def print_out_stats(fam_up, file_fam_no_graph='', wgd=''):

    """
    Prints to stdout some statistics on the families with updated synteny.

    Args:
        modified_fam (dict): a dict listing for each outgroup chromosome, the updated families
        file_fam_no_graph (str, optional): file with families that can't result in a synteny graph
        wgd (str, optional): the wgd for which the Orhtology Table was built
    """

    fam_nog = []
    if os.path.isfile(file_fam_no_graph):
        with open(file_fam_no_graph, 'r') as infile:
            fam_nog = [line.strip() for line in infile]

    fam_up_graph = [fam for fam in fam_up if fam not in fam_nog]

    print('\n')
    print("-----------------------FAMILIES WITH UPDATED SYNTENY------------------------")
    print(" Whole-genome duplication: {}".format(wgd))
    print("\n")
    print(" {} families in the final orthology table with updated synteny".format(len(fam_up)))
    if file_fam_no_graph:
        print(" For {} families, a synteny orthology graph can potentially be built"\
              .format(len(fam_up_graph)))
    print("----------------------------------------------------------------------------")
    print("\n")


def make_region_file(orthotable_file, orthotable_file_previous, corrections_file, outfile, win=15,
                     file_fam_no_graph='', wgd='', file_combin_graphs=''):
    """
    Builds and writes a file with gene families having updated synteny information in SCORPiOs
    iteration n versus iteration n-1.

    Args:
        orthotable_file (str): file with gene families at iteration n
        orthotable_file_prev (str): file with gene families at iteration n-1
        corrections_file (str): file with corrected families
        outfile (str): name of the output file.
        win (int): side of SCORPiOs sliding window for synteny orthology predictions
        file_fam_no_graph (str, optional): file with families that can't result in a synteny graph
        wgd (str, optional): the wgd for which the Orhtology Table was built
        file_combin_graphs (str, optional): summary file of graphs combination across outgroups
    """

    if not os.path.isfile(corrections_file):
        sys.stderr.write("Warning: unexisting Accepted file : {}".format(corrections_file))
        open(outfile, 'w').close()
        return

    mapping_fam_outgr = {}
    if os.path.isfile(file_combin_graphs):
        mapping_fam_outgr = read_combin_file(file_combin_graphs)

    #load genomic position of outrgoup gene for all corrected trees
    with open(corrections_file, 'r') as infile:
        corrected = [line.split('\t')[0] for line in infile]

    #load orthology table from current and previous iteration
    dcs_genes = light_load_orthotable(orthotable_file)
    dcs_genes_prev = light_load_orthotable(orthotable_file_previous)

    #create a dict of outgroup genes to keep: corrected trees or genes without ortholog in
    #SCORPiOs previous iteration
    all_modified_gene_families = get_modified_families(dcs_genes, dcs_genes_prev, corrected,
                                                       mapping_fam_outgr)


    #save all families in the same window as an updated family
    #--> fam with updated synteny information whose orthologies can be changed in next iteration
    gene_families_to_keep, genes = get_genes_to_keep(dcs_genes, all_modified_gene_families, win)


    #write regions file
    write_regions_file(gene_families_to_keep, outfile)

    #print statistics
    print_out_stats(genes, file_fam_no_graph, wgd)


def read_authorized_regions(region_file, chrom, windowsize):

    """
    Reads a file with gene families having updated synteny information in SCORPiOs iteration n
    versus iteration n-1.
    Returns a list of regions, which are bounds for windows to be considered in iteration n and a
    list of genes which are genes for which orthologies can be updated. This two list differ in
    the fact that gene can be in a considered window without having updated synteny information.

    Args:
        region_file (str): input file
        chrom (str): outgroup chromosome considered
        windowsize (int): size of the sliding window

    Returns:
        regions (list of tuple): list of regions, as tuples (start_index, stop_index),
                                 corresponding to index in the OrthologyTable.
        genes (list of str): list of genes with updated synteny information

    Note:
        If the region_file is empty, regions is set to None (and we don't filter regions in
        SCORPiOs main). If the regions file is not empty, but no family has an updated synteny
        context, regions is set to [(0, 0)] i.e. no window will be computed on this chromosome.
    """

    genes = []
    regions = []
    tmp = []

    #if region file is empty we return None
    if not os.path.getsize(region_file):
        return None, genes

    with open(region_file, 'r') as infile:

        for line in infile:

            line = line.strip().split('\t')

            #if new region
            if line == ["START"]:

                start = True

                #save current regions and start a new one
                if tmp:
                    sta = max(0, tmp[0] - windowsize)

                    #this can go out of range, will be checked in main
                    stop = tmp[1] + windowsize
                    tmp = (sta, stop)
                    regions.append(tmp)
                    tmp = []
            else:

                current_chrom, gene, ind = line

                if current_chrom == chrom:

                    #if first entry after start store start index
                    if start:

                        tmp.append(int(ind))

                    start = False

                    #update end index
                    if len(tmp) == 1:
                        tmp.append(int(ind))

                    else:
                        tmp[1] = int(ind)

                    #store gene families
                    genes.append(gene)

        #don't forget last region
        if tmp:
            sta = max(0, tmp[0] - windowsize)

            #this can go out of range, will be checked in main
            stop = tmp[1] + windowsize
            tmp = (sta, stop)
            regions.append(tmp)
            tmp = []

    #if no synteny updated family, return a region of size 0
    if not regions:
        regions = [(0, 0)]

    return regions, genes


if __name__ == '__main__':
    sys.exit()
