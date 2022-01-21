#!/usr/bin/env python

"""
    Script to build starting gene trees with TreeBeST best, from CDS back translated nucleotide
    alignments, given a species tree and a gene species mapping file.

    Example::

        $ python -m build_treebest_trees -a alis_v89.fa.gz -sp species_tree_v89.nwk
       -m genesp_v89.txt [-o treebest_forest_v89.nhx] [-nc 1] [-tmp tmp]
"""

import os
import argparse
import sys
import gzip
import multiprocessing
import traceback
import glob
import signal

from ete3 import Tree

from . import utilities as ut

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def worker_build_tree(ali, genes_sp, sptree, ali_id, tmp_folder='', X=10):

    """
    Build a gene tree from the multiple alignment string in `ali`, while accounting for the
    species tree `sptree`, using treebest best.

    If the output tree file already exists, the file will not be updated. This allows to re-
    execute a SCORPiOs snakemake run without recomputing all trees in case of error.

    Args:
        ali (str): the fasta multiple alignment
        genes_sp (str): the corresponding genes to species mapping
        sptree (str): path to the newick species tree
        ali_id (str): identifier of the tree, used in the output .nhx file name.
        tmp_folder (str): path to temp individual ali, will store temp individual tree.
        X (int, optional): -X parameter for treebest best (default=10).

    Returns:
        bool: True if no Exception was raised.

    """
    try:

        tmp_ali = tmp_folder+"tmp_ali_"+str(ali_id)+".fa"
        out_tree = tmp_folder+"tmp_tree_"+str(ali_id)+".nhx"

        if os.path.isfile(out_tree) and os.path.getsize(out_tree) > 0:
            return True

        sys.stderr.write("Building tree for alignment number "+str(ali_id)+"\n")
        sys.stderr.flush()

        mapping = {}
        genes_sp = genes_sp.strip().split('\n')
        for line in genes_sp:
            name, species = line.strip().split('\t')
            mapping[name] = species

        seq = ut.get_subali(ali, mapping, mapping)

        ut.write_fasta(seq, tmp_ali)
        cmd = "treebest best "+tmp_ali+" -f "+sptree+" -X "+str(X)+" -Z 1e-3 -q > "+out_tree
        return_value = os.system(cmd)

        #if treebest failed, we try without filtering the alignment
        #TODO it would be better to catch treebest errors using subprocess rather than os.system
        if return_value != 0:

            cmd = "treebest best "+tmp_ali+" -F 0 -f "+sptree+" -X "+str(X)+" -Z 1e-3 -q > "+\
                   out_tree
            return_value = os.system(cmd)

            if return_value != 0:
                raise Exception('treebest failed to build a tree for alignment {}'.format(ali_id))

        os.remove(tmp_ali)
        return True

    except Exception:

        traceback.print_exc()
        raise



if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-a', '--ali', type=str, help='Single file with all alignments (.fa).',
                        required=True)

    PARSER.add_argument('-sp', '--species_tree', help='Newick species tree file.', required=True)

    PARSER.add_argument('-m', '--genes_sp_map', help='Single file with corresponding gene to\
                        species mapping', required=True)

    PARSER.add_argument('-o', '--output', type=str, help='Output name for the output forest',
                        required=True)

    PARSER.add_argument('-nc', '--ncores', type=int, help='Number of threads', required=False,
                        default=1)

    PARSER.add_argument('-tmp', '--tmp_folder', type=str, help='Path for tmp invidual trees',
                        required=False, default='')

    PARSER.add_argument('-X', '--X', type=int, help='treebest -X argument', required=False,
                        default=10)

    ARGS = vars(PARSER.parse_args())

    if ARGS["tmp_folder"]:
        os.makedirs(ARGS["tmp_folder"], exist_ok=True)

    sys.stderr.write("Building starting trees with TreeBeST\n")

    OPEN = open
    if ARGS["ali"].split('.')[-1] == 'gz':
        OPEN = gzip.open

    try:
        POOL = multiprocessing.Pool(ARGS["ncores"], init_worker)

        i = 0

        with OPEN(ARGS["ali"], "rt") as INFILE_A, open(ARGS["genes_sp_map"], 'r') as INFILE_GSP:

            ASYNC_RES = []

            for i, (ALI, MAP) in enumerate(zip(ut.read_multiple_objects(INFILE_A),
                                               ut.read_multiple_objects(INFILE_GSP))):

                RES = POOL.apply_async(worker_build_tree, args=(ALI, MAP, ARGS["species_tree"], i,
                                                                ARGS["tmp_folder"], ARGS["X"]))
                ASYNC_RES += [RES]

        POOL.close()
        POOL.join()


        for RES in ASYNC_RES:
            if not RES.get():
                sys.stderr.write("An error occured in a child process\n")
                sys.exit(1)

    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        POOL.terminate()
        POOL.join()
        sys.exit(1)

    sys.stderr.write("Writing trees into a single gene tree forest file...\n")

    with open(ARGS["output"], 'w') as outforest:
        for j in range(i+1):

            TREE = Tree(ARGS["tmp_folder"]+"tmp_tree_"+str(j)+".nhx")

            #remove sp_name
            for leaf in TREE.get_leaves():
                sp_name = leaf.name.split('_')[-1]
                leaf.name = leaf.name.replace('_'+sp_name, '', 1)

            #format root node and include features and write
            TREE = TREE.write(features=["D", "S", "DD", "DCS", "B"], format_root_node=True,
                              format=1)
            outforest.write(TREE)
            outforest.write('\n//\n')

    #remove single tree files
    TMP_TREES = glob.glob(ARGS["tmp_folder"]+"tmp_tree_*.nhx")
    for tmp_tree in TMP_TREES:
        os.remove(tmp_tree)

    #remove treebest temp...
    os.remove("filtalign.fa")
