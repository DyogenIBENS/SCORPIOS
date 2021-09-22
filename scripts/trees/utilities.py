#!/usr/bin/env python

"""
    Module with functions to work with gene trees and gene alignments.

"""

import sys
import os

import numpy as np
from ete3 import Tree


def read_multiple_objects(file_object, sep="//"):

    """
    Creates a generator to read a file with several entries (trees, alignments, or other...) one
    by one.

    Args:
        file_object (file): python file object of the input file
        sep (str, optional): the separator between entries

    Yields:
        data (string) : the next tree (or alignment).
    """

    data = ""
    while True:

        line = file_object.readline()

        #don't forget to yield even if no separator at the end of file
        if not line:
            if data:
                yield data
            break

        #yield stored data each time we see the separator
        if line.strip() == sep:
            yield data
            data = ""

        else:
            data += line


def write_fasta(ali, outfile, d_names=None):

    """
    Writes a fasta alignment file.

    Args:
        ali (dict): dictionary storing the alignment {name1: seq1, name2: seq2}.
        outfile (str): name of the file to write.
        d_names (dict, optional): dictionary to transform names, for instance to add species. If
        used, all genes have to be in this dictionary.
    """

    #write to file
    with open(outfile, 'w') as outf:

        for name in ali:

            if d_names:
                assert name in d_names, "Could not find {} in genes species mapping".format(name)
                outf.write('>'+d_names[name]+'\n')

            else:
                outf.write('>'+name+'\n')

            seq = ali[name]

            #Fasta 60 characters per line
            for i in range(0, len(seq), 60):
                outf.write(seq[i:i+60]+'\n')


def delete_gaps_in_all(ali):

    """
    Removes columns of an alignment with gaps in all sequences (in-place).

    Args:
        ali (dict): dictionary storing the alignment with gene names as keys and aligned sequences
        as values.

    Note:
        Throws an assertion error if the alignment is empty

    """

    assert ali, "Empty alignment"
    lg = None
    for name in ali:
        if lg == None:
            lg = len(ali[name])
        else:
            assert len(ali[name]) == lg, "Different lengths in mutliple alignment"

    #transform string to numpy array for efficiency
    ali_array = np.array([list(ali[name]) for name in sorted(ali)])

    #define and apply a mask to remove gaps in all sequences
    mask = (ali_array == '-').all(0)
    ali_array = ali_array[:, ~mask]
    ali_list = map(''.join, ali_array)

    #transform numpy array back to dict
    for i, seq in enumerate(ali_list):
        seqname = sorted(ali)[i]
        ali[seqname] = seq


def get_subali(ali_string, genes, d_names=None):
    """
    Extract a sub-alignment of genes of interest from a fasta alignment string and
    remove columns with gaps in all sequences.

    Args:
        ali_string (str): alignment string in fasta format.
        genes (list): list of genes to extract.
        d_names (dict, optional): dictionary to add suffix to gene names, for instance to add
        species. If used, all genes have to be in this mapping dictionary.


    Returns:
        d_seq (dict): sequences sub-alignment gene names as keys and aligned sequences as values.
    """

    d_seq = {}
    extract_seq = False

    #Extract sequences present in subtree
    for line in ali_string.split('\n'):
        if line:
            if line[0] == '>':
                name = line[1:]
                if name in genes:
                    d_seq[name] = ''
                    extract_seq = True

                else:
                    extract_seq = False
            else:
                if extract_seq and '//' not in line:
                    d_seq[name] += line.strip()

    for name in genes:
        assert name in d_seq, f"Gene {name} present in tree is not found in the alignment."

    #remove gaps in all seq
    delete_gaps_in_all(d_seq)

    #add species tag
    if d_names:
        for name in list(d_seq.keys()):

            #make sure that the gene names is in the dict with new names
            assert name in d_names, "Could not find {} in genes species mapping".format(name)
            d_seq[name+'_'+d_names[name]] = d_seq.pop(name)

    return d_seq


def write_forest(inforest, outforest, corrections, current_wgd='', cor_treefiles='',
                 save_single_treefile=False):

    """
    Writes a gene tree forest after applying corrections to some gene trees, as described in
    `corrections`. Browses the input gene trees forest and writes, for each gene tree, either its
    unmodified input version or its corrected version if listed in corrections`.


    Args:
        inforest (str): Name of the .nhx file with the input gene tree forest

        outforest (str): Name of the output .nhx file for the corrected gene tree forest

        corrections (dict): For each corrected tree (described by its index in the forest)
                            a list of 5-elements tuples describing applied corrections:
                            name of wgd, corrected tree file, + 3 correction descriptors
                            If current_wgd is not used, corrections can simply be the list of
                            corrected trees, but cor_treefiles has to be specified.

        current_wgd (str, optional): If the forest is being corrected for one particular wgd, use
                                     corrections applied to this wgd. If not used, corrections
                                     apply to all trees in `corrections` but cor_treefiles has to
                                     be specified.

        cor_treefiles (str, optional): Path to corrected trees to use if current_wgd is not
                                       specified.

        save_single_treefile (bool, optional): Whether individual original input trees should be
                                               written to file and the individual corrected tree
                                               file kept.

    """

    #check args
    assert (current_wgd or cor_treefiles), "Either `current_wgd` or `cor_treefiles` should be\
                                         specified"

    stats = [0, 0, 0]
    i = -1

    with open(outforest, 'w') as outfile, open(inforest, "r") as infile:

        #for each tree in the input forest
        for i, tree in enumerate(read_multiple_objects(infile)):

            corrected = False

            #if its index is stored in correction
            if i in corrections:

                #if we correct for a particular wgd
                #check that correction stored correspond to the correct wgd
                if current_wgd:

                    wgd, cor_treefile, cor_subtrees, nb_moved_br, _ = corrections[i][-1]

                    if wgd == current_wgd:

                        corrected = True

                        stats[0] += len(cor_subtrees)
                        stats[1] += 1
                        stats[2] += len([i for i in nb_moved_br if i])

                #if we correct for all use the cor_treefiles arg to locate the corrected tree file
                else:
                    corrected = True
                    cor_treefile = cor_treefiles+str(i)+'.nhx'

                #write corrected tree
                if corrected:
                    with open(cor_treefile, 'r') as treef:
                        for line in treef:
                            outfile.write(line)

                    #write input tree if requested
                    if save_single_treefile:
                        wtree = Tree(tree, format=1)
                        cor_folder = '/'.join(cor_treefiles.split('/')[:-1])
                        wtree.write(outfile=cor_folder+'/ori_'+str(i)+'.nhx', format=1,
                                    features=["S", "D", "DD", "DCS"], format_root_node=True)

                    #if not requested remove also corrected version
                    else:
                        os.remove(cor_treefile)

            if not corrected:
                outfile.write(tree)

            outfile.write('\n//\n')

    if current_wgd:
        print('\n')
        print("----------------------------SUBTREES RE-GRAFTING----------------------------")

        print(" Whole-genome duplication: {}".format(current_wgd))

        print('\n')

        print(" On {} total gene trees, {} total corrected subtrees in {} different trees"\
              .format(i+1, stats[0], stats[1]))

        print(" {} corrections required topological changes for other branches in the tree"\
              .format(stats[2]))

        print("----------------------------------------------------------------------------")
        print('\n')


if __name__ == '__main__':
    sys.exit()
