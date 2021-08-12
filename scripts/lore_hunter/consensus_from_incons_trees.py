#! /usr/bin/env python

import glob

import dendropy
from ete3 import Tree

def filter_trees_strict(input_trees, outgr):
    to_keep = []
    for input_tree in input_trees:
        tree = Tree(input_tree)
        all_sp = {i.S for i in tree.get_leaves()}.difference(outgr)

trees = dendropy.TreeList()



for tree_file in ['pythonidae.mb.run1.t',
        'pythonidae.mb.run2.t',
        'pythonidae.mb.run3.t',
        'pythonidae.mb.run4.t']:
    trees.read_from_path(
            tree_file,
            'nexus',
            tree_offset=20)
con_tree = trees.consensus(min_freq=0.5)
print(con_tree.as_string('newick'))