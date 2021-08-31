#!/usr/bin/env python

"""
prototype, quick and dirty-made script to improve if it proves useful
"""
import sys

import glob
import itertools
import os
import multiprocessing

import argparse
import random
import numpy as np

from ete3 import Tree

import tree_distance as td

#FIXME: error handling in multiprocessing

def load_incons(finput):
    with open(finput, 'r') as infile:
        incons = []
        incons_all = []
        for line in infile:
            _, genes, cat = line.strip().split('\t')
            if "Incons" in cat:
                incons_all += genes.split()
            if cat == "Inconsistent":
                incons += genes.split()
    return set(incons), set(incons_all)

def product(*iterables, **kwargs):
    if len(iterables) == 0:
        yield ()
    else:
        iterables = iterables * kwargs.get('repeat', 1)
        myiter = iterables[0]
        for item in myiter() if callable(myiter) else iter(myiter):
            for items in product(*iterables[1:]):
                yield (item, ) + items


def all_trees_to_newick(tree_list):
    trees_dict = {}

    for tree in tree_list:
        #FIXME do smthg better for this
        treename = os.path.splitext(os.path.basename(tree))[0].replace("_final", "")
        mytree = Tree(tree)
        leaves = mytree.get_leaves()
        species = {i.S.split("_")[0] for i in leaves}

        if len(species) > 2:

            for leaf in leaves:

                leaf.name = leaf.S.replace('_', '')

            trees_dict[treename] = mytree

    return trees_dict


def apply_perm(tree, perm, leaves):

    perm = list(perm)

    d_convert = {}

    for group in perm:
        k = 0
        for leaf in leaves:
            species = leaf[:-1]

            if species in [i[:-1] for i in group]:

                d_convert[leaf] = group[k]
                k += 1

    for leaf in tree.get_leaves():

        leaf.name = d_convert[leaf.name]

    return tree

def get_all_permutation(tree):

    dsp = {}
    leaves_names = [i.name for i in tree.get_leaves()]
    for leaf in leaves_names:
        species = leaf[:-1]
        dsp[species] = dsp.get(species, [])
        dsp[species].append(leaf)

    all_permutations = []
    for leaves in dsp.values():
        all_permutations.append(list(itertools.permutations(leaves)))

    return itertools.product(*all_permutations), leaves_names

def compare_trees(pair, res, dist_method="Euclid"):

    tree1, tree2 = pair
    tree1_name, tree1_ete3 = tree1
    tree2_name, tree2_ete3 = tree2


    mytree1 = tree1_ete3.copy()
    mytree2 = tree2_ete3.copy()

    min_dist = np.inf

    leaves_1 = {i.name for i in mytree1}
    leaves_2 = {i.name for i in mytree2}

    intersect = leaves_1.intersection(leaves_2)

    if len(leaves_2) < len(leaves_1):
        mytree1, mytree2 = mytree2.copy(), mytree1.copy()
        tree1_name, tree2_name = tree2_name, tree1_name

    all_permut_t1, freeze_leaves1 = get_all_permutation(mytree1)


    if len(intersect) > 2:

        for perm in all_permut_t1:

            t1_perm = apply_perm(mytree1.copy(), perm, freeze_leaves1)
            t2_perm = mytree2.copy()

            leaves_1 = {i.name for i in t1_perm}
            leaves_2 = {i.name for i in t2_perm}

            intersect = leaves_1.intersection(leaves_2)

            if len(intersect) > 2:

                t1_perm.prune(intersect, preserve_branch_length=True)
                t2_perm.prune(intersect, preserve_branch_length=True)

                t1_nwk = t1_perm.write(format=1)
                t2_nwk = t2_perm.write(format=1)


                if dist_method == "Euclid":
                    dist = td.getEuclideanDistance(td.PhyloTree(t1_nwk, True),\
                                                   td.PhyloTree(t2_nwk, True), True)

                elif dist_method == "RF":

                    dist = td.getRobinsonFouldsDistance(td.PhyloTree(t1_nwk, True),\
                                                        td.PhyloTree(t2_nwk, True), True)

                else:
                    sys.stderr.write("Error: unsupported dist_method '"+dist_method+\
                                     "' accepted values are 'Euclid' and 'RF'\n")
                    sys.exit(1)


                if dist < min_dist:

                    min_dist = dist

    if  min_dist == np.inf:
        min_dist = np.nan

    res[(tree1_name, tree2_name)] = min_dist


def grouper(size, iterable, fillvalue=None):
    "grouper(size, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * size
    return itertools.izip_longest(fillvalue=fillvalue, *args)


if __name__ == '__main__':


    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--treesdir', required=True)

    PARSER.add_argument('-i', '--incons', required=True)

    PARSER.add_argument('-o', '--outdir', required=True)

    PARSER.add_argument('-nc', '--ncores', required=False, default=1, type=int)

    PARSER.add_argument('-r', '--ratio', required=False, default=1.3, type=int)

    PARSER.add_argument('-m', '--max_incons_trees', required=False, default=2000, type=int)

    ARGS = vars(PARSER.parse_args())

    trees = sorted(glob.glob(ARGS["treesdir"]+'/*_final.nhx'))

    ALL_TREES = all_trees_to_newick(trees)

    print(len(ALL_TREES))

    INCONS, IALL = load_incons(ARGS["incons"])
    print(len(INCONS), len(IALL))

    for t in ALL_TREES:
        leaves = {i.name for i in ALL_TREES[t].get_leaves()}
        intersect = INCONS.intersection(leaves)
        if intersect:
            new_key = list(intersect)[0]
            ALL_TREES[new_key] = ALL_TREES.pop(t)

    incons_trees = [i for i in ALL_TREES if i in INCONS]
    other_trees = [i for i in ALL_TREES if i not in IALL]
    print(len(incons_trees), len(other_trees))

    #FIXME Ratio to give as a param
    nb_incons = min(int(0.95*len(incons_trees)), ARGS["max_incons_trees"])
    SAMPLE_INCONS = random.sample(incons_trees, ))

    #random.sample(other_trees, int(5000-0.9*len(incons_trees)))
    tmp = int(0.95*len(incons_trees))
    size_others = min(int(tmp*ARGS["ratio"]), len(other_trees))
    SAMPLE_OTHER = random.sample(other_trees, size_others)

    # SAMPLE_INCONS = incons_trees
    # SAMPLE_OTHER = other_trees

    print(len(SAMPLE_OTHER), len(SAMPLE_INCONS))
    ALL_TREES = {k:ALL_TREES[k] for k in SAMPLE_INCONS+SAMPLE_OTHER}
    trees_pairs = itertools.combinations(ALL_TREES.iteritems(), 2)

    loop = 1

    print(len(ALL_TREES))
    TOT = len(ALL_TREES)*(len(ALL_TREES)-1)/2

    GROUP_BY = 10000

    OUTDIR = ARGS["outdir"]

    for tree_pair_subset in grouper(GROUP_BY, trees_pairs):

        MANAGER = multiprocessing.Manager()
        RES = MANAGER.dict()
        POOL = multiprocessing.Pool(ARGS["ncores"])
        # ASYNC_RESULTS = []
        for pair in tree_pair_subset:

            if pair:

                POOL.apply_async(compare_trees, args=(pair, RES))

        POOL.close()
        POOL.join()

        sys.stderr.write("Computed "+str(max(loop * GROUP_BY, TOT)) +" distances out of "+str(TOT)+\
                         " ("+str(max(round(loop*GROUP_BY/float(TOT) * 100, 2)), 100)+" %) \n")

        with open(OUTDIR+"/dist_mat_"+str(loop)+".csv", 'w') as fw:
            for comp in RES.keys():

                min_dist = RES[comp]

                fw.write(comp[0]+','+comp[1]+','+str(min_dist)+'\n')

        loop += 1
