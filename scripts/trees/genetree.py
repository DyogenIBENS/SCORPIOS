#!/usr/bin/env python

"""
    Module with functions to work with a gene tree.
"""

import sys

from itertools import chain

import numpy as np
from ete3 import Tree


def load_corrections(files):

    """
    Gets the name and path of SCORPiOs corrected subtrees (i.e accepted by AU-tests).

    Arg:
        files (str): Comma-delimited list of files with accepted corrections.

    Returns:
        dict: For each corrected subtree, the path to the tree file and the name of the corrected
              WGD.
    """

    cor_subtrees = {}

    for input_file in files.split(','):

        with open(input_file, 'r') as infile:

            for line in infile:

                line = line.strip().split('\t')
                outgr_gene, tree_path, wgd = line[0:3]
                cor_subtrees[outgr_gene] = (tree_path, wgd)

    return cor_subtrees


def get_solution_subtree(corr_dict, subtree_name):

    """
    Loads a SCORPiOs-corrected subtree.

    Args:
        corr_dict (dict): SCORPiOs-corrected subtrees
        subtree_name (str): name of the corrected subtree to load

    Returns:
        ete3.Tree: gene tree in an ete3.Tree object
    """

    outgr_gene = '_'.join(subtree_name.split('_')[:-1])

    soldir, wgd = corr_dict[outgr_gene]

    tree_name = soldir+'/'+wgd+'/'+outgr_gene+'.nh'

    sys.stderr.write('Correction using solution subtree: {}\n'.format(tree_name))

    stree = Tree(tree_name)

    stree.set_outgroup(stree&subtree_name)

    return stree


def save_nhx_tags(attribute_names, groups_of_leaves):

    """
    Stores node features in a dictionary. This dictionary can be used to add .nhx
    tags to a tree.

    Args:
        attribute_names (list of str): list of names of features to add
        groups_of_leaves (nested list): For each feature, sub-lists of nodes to annotate with the
                                        same feature value, the value being the index of the
                                        sub-list.
    Returns:
        dict: For each node, the attribute name and its associated value.

    Example:
        if `attribute_names` is ['a', 'b'] and `groups_of_leaves` is
        [[['gene1', 'gene2'], ['gene3']], ['gene4']], the function returns:
        features = { 'gene1': [('a', 1)], 'gene2': [('a', 1)], 'gene3': [('a', 2)],
                     'gene4': [('b', 1)] }
    """

    features = {}

    for i, groups in enumerate(groups_of_leaves):

        for j, group in enumerate(groups):

            for node in group:

                features[node] = features.get(node, [])

                features[node].append((attribute_names[i], str(j+1)))

    return features


def add_nhx_tags(wtree, features):

    """
    Adds leaves attribute stored in a dict as .nhx tags (i.e as features in the ete3.Tree object).
    Leaf names in dict can be with or without species tags, but leaves should have name and S
    attributes.

    Args:
        wtree (ete3.Tree): Tree for which to add new leaves attributes (modified in-place)
        features (dict): A dictionary giving, for each leaves, the pairs (attribute, value) to add.

    Returns:
        list: the list of names of all added features
    """

    added_features = []
    for leaf in wtree.get_leaves():

        gene_name_sp = leaf.name+'_'+leaf.S
        leaf_features = []

        if leaf.name in features:
            leaf_features = features[leaf.name]

        elif gene_name_sp in features:
            leaf_features = features[gene_name_sp]


        if leaf_features:
            for leaf_feature in leaf_features:
                attribute, value = leaf_feature
                leaf.add_feature(attribute, value)

                if attribute not in added_features:
                    added_features.append(attribute)
    return added_features



def branch_length_closest(tree, gene, group_of_genes):

    """
    Finds the gene closest to `gene` in `tree` amongst a `group_of_genes`, i.e with the shortest
    branch length.

    Args:
        tree (ete3.Tree): input tree
        gene (ete3.TreeNode): gene for which to search a neighbour
        group_of_genes (list of str): list of candidate neighbour genes

    Returns:
        (str): name of the closest neighbour, in terms of branch lengths
    """

    dist_min = np.inf
    min_gene_d = {}
    names = {}

    #iterate over genes in `group_of_genes` and compute branch-length distance
    for target in group_of_genes:

        dist = tree.get_distance(gene, target)

        if dist <= dist_min:
            dist_min = dist
            min_gene_d[target.name] = dist_min
            names[target.name] = target

    #arbitrary choice to ensure deterministic answer
    all_max_genes = []
    for hit_gene in min_gene_d:
        if min_gene_d[hit_gene] == dist_min:
            all_max_genes.append(hit_gene)
    best_gene = sorted(all_max_genes)[0]
    return names[best_gene]


def closest_gene_in_tree(tree, node, group_of_genes, attr='name'):

    """
    Finds the gene closest to `gene` in `tree` amongst a `group_of_genes`, i.e the gene with
    (i) the shortest topological distance to `gene` and (ii) the shortest branch-length distance
    in case of ties in (i).

    Args:
        tree (ete3.Tree): input tree
        gene (ete3.TreeNode): gene for which to search for closest neighbour
        group_of_genes (list of str): list of candidate neighbour genes
        attr (str, optional): name of the attriute storing gene names


    Returns:
        (str): name of closest neighbour
    """

    sp_gene = ''

    #we go up the tree starting from `node`
    while node.up:

        node = node.up

        #genes at this topological distance
        descendant_genes = [i for i in node.get_leaves() if getattr(i, attr) in group_of_genes]

        if descendant_genes:

            #if no tie in topological distance, we found closest neighbour
            if len(descendant_genes) == 1:
                sp_gene = descendant_genes[0]

            #else use branch-lengths
            else:
                sp_gene = branch_length_closest(tree, node, descendant_genes)

            break


    return sp_gene


def find_node_with_most_desc(tree1, corrected_leaves):

    """
    Finds the node in `tree1` with all its descending leaves in `corrected_leaves` and the maximum
    number of descending leaves.

    Args:
        tree1 (ete3.Tree): input tree
        corrected_leaves (list of str): list of leaves name to maximize in descendants

    Returns:
        ete3.TreeNode: the identified node

    """

    node2leaves = tree1.get_cached_content(store_attr='name') #we store a light version of the tree
    inter_max = 0
    for node in tree1.traverse():

        descendants = set(node2leaves[node])

        #if all descending leaves are in corrected_leaves
        if not descendants.difference(corrected_leaves):

            #we maximize the number of descendants
            inter = len(descendants)
            if inter > inter_max:
                node_max = node
                inter_max = inter

    return node_max


def find_sister_of_outgroup(leaf_outgr, authorized_sp, sister_outgroup_genes):

    """
    Extracts a list of genes related to an outgroup gene `leaf_outgr`: genes belonging to related
    species `authorized_sp` and grouped together in the gene tree.

    Args:
        leaf_outgr (ete3 TreeNode): node of the outgroup gene in the tree
        authorized_sp (list of str): list of related species
        sister_outgroup_genes (list of str): list of related genes (to update in-place)
    """

    is_sister = True

    #we go up the tree from outgr_leaf
    while is_sister and leaf_outgr.up:

        leaf_outgr = leaf_outgr.up

        #genes below the node
        new_gene = [i for i in leaf_outgr.get_leaves() if i.name not in sister_outgroup_genes]

        #we continue while all genes are in `authorized_sp`
        if new_gene:
            if new_gene[0].S not in authorized_sp:
                is_sister = False


        #store related genes in-place
        if new_gene:
            for gene in new_gene:
                if gene.S in authorized_sp:
                    sister_outgroup_genes.append(gene.name)


def keep_sis_genes_together(duplicated_sp_subtree, outgr, sister_outgroup_genes, outgroup_subtree,
                            node_max='node_max'):
    """
    Keeps genes of all outgroup species together when modifying a gene tree, so that the
    new tree remains species tree consistent for these species that branch between the outgroup and
    duplicated species.

    Args:
        duplicated_sp_subtree (ete3.Tree) : Synteny-corrected subtree
        outgr (str) : name of the non-duplicated outgroup gene
        sister_outgroup_genes (list of str) : genes that are grouped with the outgroup gene in the
                                              original tree and in related species
        outgroup_subtree (ete3.Tree) : subtree with only outgroup and related genes
        node_max (str, optional) : internal node name in the outgroup subtree where to paste the
                                   duplicated species subtree. If empty a new tree combining both
                                   is created.

    Returns:
        ete3.Tree : a new tree where the outgroup gene in the synteny-corrected is replaced by the
                    subtree of all outgroup genes
    """

    outgr_subtree_leaves = [i.name for i in outgroup_subtree.get_leaves() if i.name in\
                                                                             sister_outgroup_genes]
    leaves_to_prune = [i.name for i in duplicated_sp_subtree.get_leaves() if i.name != outgr]

    if node_max:
        outgroup_subtree.prune(outgr_subtree_leaves + [node_max])
        node = outgroup_subtree.search_nodes(name=node_max)[0]
        node.name = ''
        duplicated_sp_subtree.prune(leaves_to_prune)
        node.add_child(duplicated_sp_subtree.copy())
        duplicated_sp_subtree = outgroup_subtree.copy()

    else:
        outgroup_subtree.prune(outgr_subtree_leaves)
        duplicated_sp_subtree.prune(leaves_to_prune)
        new = Tree()
        new.add_child(outgroup_subtree)
        new.add_child(duplicated_sp_subtree.copy())
        duplicated_sp_subtree = new

    return duplicated_sp_subtree


def keep_subsequent_wgd_species(stree, ensembl_tree, missing_leaves_keep, sp_current_wgd,
                                authorized_sp):

    """
        When re-grafting a subtree corrected for species descending from 'WGD1' only, keep
        positions of species with subsequent WGDs consistent in the tree.
        To do so, find the closest 'WGD1-only' species gene in ensembl tree and keep subsequently
        duplicated species genes at the same position (relative to it).

        Modifies `stree` in-place.

        Args:
            stree (ete3.Tree): Tree object for the synteny corrected tree of WGD1
            ensembl_tree (ete3.Tree): Tree object for the full original gene tree
            missing_leaves_keep (list of ete3.TreeNode): Genes of subsequently duplicated species
            sp_current_wgd (list of str): List of WGD1 duplicated species
            authorized_sp (dict): Dict used to keep the tree consistent with the species tree. For
            a 'WGD1' species, a list of WGD1 species that are closer to it than are 4R species.
    """

    #genes in the WGD1 corrected tree
    sleaves = [i.name for i in stree.get_leaves()]
    stree.prune([i for i in stree.get_leaves()])

    #genes in subsequent WGDs
    missing = [i.name for i in missing_leaves_keep]

    #original tree
    enstree = ensembl_tree.copy()

    #find all clades of subsequent WGD genes to replace them at a correct position together
    for leaf in enstree.get_leaves():
        if leaf.name in missing:
            leaf.missing = "Y"
    clades = enstree.get_monophyletic(values=["Y"], target_attr="missing")

    closest_gene = {}

    #for each clade to replace
    for clade in clades:

        #find closest neighbour in the original tree which is in WGD1 duplicated species
        outgr_gene = closest_gene_in_tree(enstree, clade, sp_current_wgd, attr='S')

        if outgr_gene.name in sleaves:
            closest_gene[clade] = outgr_gene

    #if the closest WGD1 gene is in the WGD1 stree
    #we'll keep 4R genes in the family, at a similar position
    for outgr_gene in set(closest_gene.values()):
        clades = [i for i in closest_gene if closest_gene[i] == outgr_gene]
        clades = list(chain.from_iterable(clades))

        subtree_4r = enstree.copy()
        subtree_4r.prune([i.name for i in clades]+[outgr_gene.name])
        outgroup_subtree = stree.copy()

        outgr = stree.get_leaves_by_name(name=outgr_gene.name)[0]

        sister_outgroup_genes = [outgr_gene.name]
        find_sister_of_outgroup(outgr, authorized_sp[outgr_gene.S], sister_outgroup_genes)

        #We keep all sister outgroup genes together in the corrected tree
        if len(sister_outgroup_genes) > 1:

            #stree is modified in-place
            subtree_4r = keep_sis_genes_together(subtree_4r, outgr_gene.name,
                                                 sister_outgroup_genes, outgroup_subtree,
                                                 node_max='')
            lca = stree.get_common_ancestor(sister_outgroup_genes)
            # cop = stree.copy()
            stree.prune([lca] + [i for i in stree.get_leaves()\
                                 if i.name not in sister_outgroup_genes])
        else:
            lca = outgr

        #in case we do not paste the subtree at a terminal node (cleared above)
        if len(lca.children) == 2:
            lca_cop = lca.copy()
            tmp = Tree()
            lca_cop.prune([i for i in lca_cop.get_leaves()])
            tmp.add_child(lca_cop.copy())
            tmp.add_child(subtree_4r.copy())
            lca.up.add_child(name='here')
            lca.detach()
            lca_new = stree.search_nodes(name="here")[0]
            lca_new.add_child(tmp.copy())
            lca_new.name = ''

        else:
            lca.name = ''
            lca.add_child(subtree_4r)

    #remove potential signle child internal nodes artefact
    stree.prune([i for i in stree.get_leaves()])

    #clean up attributes
    for leaf in stree.get_leaves():
        if hasattr(leaf, 'missing'):
            delattr(leaf, 'missing')


def copy_nhx_tags(tree_ref_tags, tree_target):

    """
    Copies nhx tags stored in leaves of tree1 to leaves of tree2. tree2 is modified in-place

    Args:
        tree_ref_tags (ete3.Tree) : tree with nhx tags to copy
        tree_target (ete3.Tree) : tree to copy tags to

    """

    d_tags = {}
    for leaf in tree_ref_tags.get_leaves():
        attr = {i:getattr(leaf, i) for i in vars(leaf) if i not in\
                                   ["dist", "name", "support", "features"] and i[0] != '_'}
        d_tags[leaf.name] = attr

    for leaf in tree_target.get_leaves():
        if leaf.name in d_tags:
            for tag in d_tags[leaf.name]:
                value = d_tags[leaf.name][tag]
                setattr(leaf, tag, value)



if __name__ == '__main__':
    sys.exit()
