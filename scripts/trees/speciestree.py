#!/usr/bin/env python

"""
    Module with functions to work with a species tree.
"""

import sys
from collections import OrderedDict
from ete3 import Tree


def search_one_node(tree, node_name):

    """
    Searches for a node in the input tree given its name. Throws AssertionError if the node name
    is not found or is not unique.

    Args:
        tree (ete3 Tree): input tree
        node_name (str): node to search

    Returns:
        ete3 TreeNode: the matched node
    """

    nodes = tree.search_nodes(name=node_name)
    assert nodes, "{} not in tree".format(node_name)
    assert len(nodes) == 1, "{} ambiguous node name".format(node_name)
    return nodes[0]


def remove_anc(tree_file, out_file):

    """
    Removes any internal node name, such as ancestor names, in the input tree and writes it to a
    new file.

    Args:
        tree_file (str): Path to the input newick formatted tree.
        out_file (str): Path for the output file.
    """

    tree = Tree(tree_file, format=1)
    tree.prune([i for i in tree.get_leaves()])
    tree.write(outfile=out_file, format=9)


def get_species(species_tree, anc, other_wgd_anc='', lowcov_species=''):

    """
    Extracts a list of species descending from a given ancestor in a species tree. Filter out
    species under particular ancestors (i.e subsequent WGDs for instance) given by
    'other_WGD_anc', as well as 'low-coverage' species given by 'lowcov_species'.

    Args:
        tree_file (str): Path to the input newick tree.
        out_file (str): Path for the output file.
        other_wgd_anc (str, optional): Comma-delimited names of ancestors with subsequent WGDs.
        lowcov_species (str, optional): Comma-delimited names of 'lowcoverage' species to exclude.

    Returns:
        species (list): The list of species.
    """

    #get species under the WGD node
    tree = Tree(species_tree, format=1)
    lca = search_one_node(tree, anc)

    species = [i.name for i in lca.get_leaves()]

    #get species under subsequent WGD nodes
    if other_wgd_anc:
        other_wgd_anc = other_wgd_anc.split(',')
        if anc not in other_wgd_anc:
            other_wgd_anc.append(anc)
        other_wgd_anc = get_anc_order(species_tree, other_wgd_anc, tips_to_root=False)
        all_wgd = other_wgd_anc[anc] #WGD that occurred after 'anc'
        species_with_other_wgd = []
        for wgd in all_wgd:
            lca = search_one_node(tree, wgd)
            species_with_other_wgd.extend([i.name for i in lca.get_leaves()])

        species = [sp for sp in species if sp not in species_with_other_wgd]

    #filter out 'lowcov'
    if lowcov_species:
        list_lowcov_sp = lowcov_species.split(',')
        species = [sp for sp in species if sp not in list_lowcov_sp]

    return species


def get_sister_species(species_tree, species, anc):

    """
    Extracts a list of species related to a given species: species branching between `species` and
    the ancestor `anc`.

    Args:
        species_tree (ete3 Tree): ete3 tree object
        species (str): name of the species
        anc (str): name of the ancestor

    Returns:
        list: species branching between `species` and `anc`
    """

    sp_and_sisters = [species]
    duplicated_sp = get_species(species_tree, anc)
    tree = Tree(species_tree, format=1)
    lca = tree.get_common_ancestor([species]+duplicated_sp)
    sp_and_sisters += [i.name for i in lca.get_leaves() if i.name not in [species]+duplicated_sp]

    return sp_and_sisters


def is_below(node1, node2):

    """
    Checks if node2 is below node1 in the tree topology.

    Args:
        node1 (ete3 TreeNode): node1
        node2 (ete3 TreeNode): node2

    Returns:
        bool: True if node2 is below node1, False otherwise.
    """

    below = False
    all_below = [i.name for i in node1.get_descendants()]
    if node2 in all_below:
        below = True
    return below


def get_anc_order(tree_file, ancestors=None, tips_to_root=False, prune=True):

    """
    Orders input ancestors with respect to their position in the species tree. Can be ordered from
    root to tips (default) or tips to root.

    Args:
        tree_file (str): Path to the input newick formatted tree.
        ancestors (optional, list of str): List of ancestor names. If unspecified, all the ancestors
                   in the trees will be returned.

    Returns:
        OrderedDict: ancestor names in the requested order (keys) and list of ancestors in the
        input list that are below it (values).
    """

    tree = Tree(tree_file, format=1)
    if not ancestors:
        ancestors = [i.name for i in tree.traverse() if not i.is_leaf()]
    if prune:
        tree.prune([i for i in tree.get_leaves()])
    dist_to_root = {i:tree.get_distance(i) for i in ancestors}
    anc_order = sorted(dist_to_root, key=dist_to_root.get)

    if tips_to_root:
        anc_order = anc_order[::-1]

    anc_order_dict = OrderedDict()
    for anc in anc_order:

        anc_order_dict[anc] = []
        anc_node = search_one_node(tree, anc)

        for anc2 in ancestors:

            if anc != anc2:
                if is_below(anc_node, anc2):
                    anc_order_dict[anc].append(anc2)

    return anc_order_dict


if __name__ == '__main__':
    sys.exit()
