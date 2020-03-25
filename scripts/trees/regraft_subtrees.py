#!/usr/bin/env python

"""
    Script to re-graft corrected subtrees in their original tree and write the corrected
    gene trees forest.

    Example:
        $ python -m scripts.trees.regraft_subtrees -t trees.nhx -a alis.fa -s species_tree.nwk
                                                   -acc Accepted_trees.txt -o outtrees.nhx
                                                   -anc Clupeocephala,Salmonidae
                                                   -ogr Lepisosteus.oculatus,Amia.calva_Esox.lucius
                                                   [-n 1] [-tmp path/tmp] [-sa n] [-br y]
"""

import multiprocessing
import traceback
import os
import sys
import argparse
import gzip
import itertools

from ete3 import Tree


from . import utilities as ut
from . import genetree as gt
from . import speciestree as spt


def topo_changes(lca, stree, leaves_to_move, outgr, authorized_sp):

    """
    Makes necessary topological changes in the initial tree to paste a corrected subtree. Genes
    in the corrected subtree are grouped as a clade in the new tree, while other branches are
    modified as less as possible.

    Args:
        lca (ete3.TreeNode): the original tree topology below the node containing all corrected
                             leaves
        stree (ete3.Tree): the corrected subtree to re-graft
        leaves_to_move (list of str): genes that split genes of the corrected subtrees into
                                      several clades the original tree
        outgr (str): Name of the outgroup gene used to build the corrected subtree
        authorized_sp (list of str): list of species whose genes should remained grouped with the
                                     outgroup

    Returns:
        ete3.Tree: The new gene tree with modified topology
    """

    #the final tree will consist in three subtrees pasted together:
    #the corrected duplcated_species_subtree, outgroup_genes and genes to move outside (final_tree)
    duplicated_sp_subtree = stree.copy("newick-extended")
    final_tree = lca.copy("newick-extended")

    #find node with max corrected trees descendants
    #node_max is the position where to paste the corrected subtree into the original tree
    leaves_in_cor = {i.name for i in stree.get_leaves()}
    node_max = gt.find_node_with_most_desc(final_tree, leaves_in_cor)
    if not node_max.is_leaf():
        node_max.name = 'node_max'

    outgroup_subtree = final_tree.copy("newick-extended")

    #find all sister genes of outgr that are together in the tree
    #(we'll keep them together in the tree)
    leaf_outgr = [i for i in final_tree.get_leaves() if i.name == outgr][0]
    sister_outgroup_genes = [leaf_outgr.name]
    gt.find_sister_of_outgroup(leaf_outgr, authorized_sp, sister_outgroup_genes)
    sister_outgroup_genes = [i for i in sister_outgroup_genes if i not in\
                            leaves_in_cor or i == leaf_outgr.name]

    #genes not in the constraint and not outgroup genes will be kept like they are but as outside
    #the tree topology for these genes is stored in final_tree
    pruned = False
    non_outgroup_leaves = [i for i in leaves_to_move if i not in sister_outgroup_genes]
    if non_outgroup_leaves:
        pruned = True
        final_tree.prune([i for i in non_outgroup_leaves] + [node_max])

    #We keep all sister outgroup genes together and as direct outgroups of the corrected tree
    if len(sister_outgroup_genes) > 1:
        #duplicated_sp_tree is modified in-place
        duplicated_sp_subtree = gt.keep_sis_genes_together(duplicated_sp_subtree, outgr,
                                                           sister_outgroup_genes,
                                                           outgroup_subtree, node_max.name)

    #if we have remaining leaves (genes outside)
    #we paste the duplicated_sp_subtree+outgroups at node_max
    if pruned:
        node_max.add_child(duplicated_sp_subtree)
        node_max.name = ''

    #otherwise we have placed all leaves as direct outgr and we don't need to update
    else:
        final_tree = duplicated_sp_subtree.copy("newick-extended")

    #Finally, we remove artefactual single-child internal nodes (due to the mutliple copy-pasting)
    final_tree.prune([i.name for i in final_tree.get_leaves()])

    return final_tree, sister_outgroup_genes


def correct_wtrees(tree, to_cor, res, tree_id, outfiles, outgroup_sp, sp_below_wgd=None,
                   sp_current_wgd=None, tag=''):

    """
    Re-graft all subtrees corrected for one WGD into the initial gene tree.

    Args:
        tree (ete3.Tree): the initial gene tree
        to_cor (dict): for each corrected subtree (identified by the outgroup gene), the path to
                       the corrected subtree (value)
        res (dict): Stores a summary of applied corrections for `tree`
        tree_id (int): Index of the gene tree in the forest, used as key in `res`
        outfiles (str): path to store the gene tree after subtree re-grafting
        outgroup_sp (dict): for each outgroup and WGD species, a list of
                            its sister species
        sp_below_wgd (dict, optional): For each subsequent wgds, the duplicated species and their
                                       outgroups
        sp_current_wgd (list, optional): A list of duplicated species for the WGD for which the
                                         subtrees were corrected
        tag (str, optional): tag to add (WGD for instance) to document corrections in outputs

    Note:
        The `res` dictionary is filled in-place with a correction summary for `tree`.
    """

    wtree = Tree(tree, format=1)

    #all outgroup genes in the tree
    subtrees = [i for i in wtree.get_leaves() if i.S in outgroup_sp.keys()]

    #all corrected subtrees to re-graft
    cor_subtrees = [i for i in subtrees if i.name in to_cor]
    if not cor_subtrees:
        return

    else:
        whole_tree = ('_').join([i.name for i in cor_subtrees])

        d_sp = {}
        for leaf in wtree.get_leaves():
            leaf.name = leaf.name+'_'+leaf.S
            d_sp[leaf.name] = leaf.S

        all_missing_leaves, all_subtrees_leaves = [], []

        #Re-graft each subtree that we corrected
        for cor_subtree in cor_subtrees:
            missing_leaves = []
            corrected_tree = gt.get_solution_subtree(to_cor, cor_subtree.name)
            cor_descendants = []
            for leaf in corrected_tree.get_leaves():
                leaf.add_features(S=d_sp[leaf.name])
                cor_descendants.append(leaf.name)

            all_subtrees_leaves.append(cor_descendants)
            lca = wtree.get_common_ancestor(cor_descendants)
            leaves_under_lca = [i.name for i in lca.get_leaves()]

            #need to modify the topology of the tree if some genes from other species
            #are present under the node to correct and split the family
            if set(leaves_under_lca) != set(cor_descendants):

                #check for leaves under the node that are not in the corrected tree
                missing_leaves = [i for i in leaves_under_lca if i not in cor_descendants]

                #keep species below subsequent WGDs at a consistent place
                if sp_below_wgd:
                    for subs_wgd in sp_below_wgd:
                        leaves_to_keep_below = [i for i in lca.get_leaves()\
                                        if i.name not in cor_descendants\
                                        and i.S in sp_below_wgd[subs_wgd]]
                        gt.keep_subsequent_wgd_species(corrected_tree, lca, leaves_to_keep_below,
                                                       sp_current_wgd, outgroup_sp[subs_wgd])

                cor_descendants = [i.name for i in corrected_tree.get_leaves()]
                missing_leaves = [i for i in leaves_under_lca if i not in cor_descendants]

                #make the topological changes necessary to place the subtree
                if missing_leaves:
                    corrected_tree, outgroup_leaves = topo_changes(lca, corrected_tree,\
                                                                   missing_leaves,\
                                                                   cor_subtree.name,\
                                                                   outgroup_sp[cor_subtree.S])

            all_missing_leaves.append([i for i in missing_leaves if i not in outgroup_leaves])

            #replace node if it is not root
            if lca.up:
                node = lca.up
                lca.detach()
                node.add_child(corrected_tree)

            #otherwise directly replace tree
            else:
                wtree = corrected_tree.copy("newick-extended")

        #save edition tags in a dict, we will put them as NHX attributes in the end
        #otherwise they would be wiped out by treebest during reconciliation or br-length computing
        #This should be checked again, I am not sure treebest removes .nhx attributes anymore
        #Perhaps storing attributes and putting them back at the end is not nescessary
        d_edit = gt.save_nhx_tags(['CORR_ID_'+tag, 'MOVED_ID_'+tag],
                                  [all_subtrees_leaves, all_missing_leaves])

        for leaf in wtree.get_leaves():
            leaf.name = leaf.name.replace('_'+leaf.S, '', 1)

        final_wtree_file = outfiles+whole_tree
        wtree.write(outfile=final_wtree_file, format=1,
                    features=["S"], format_root_node=True)
        size_of_recalc_trees = [len(i) for i in all_missing_leaves]
        res[tree_id] = res.get(tree_id, [])
        res[tree_id].append((tag, final_wtree_file, cor_subtrees, size_of_recalc_trees, d_edit))




def worker_rec_brlgth(tree, outfolder, treeid, sptree, ali='', prefix='cor', corrections=None,
                      brlengths=True, resume=False):

    """
    Reconciles a given gene tree with the species tree using treebest sdi, and optionaly computes
    branch-length using treebest phyml. Also adds .nhx tags to corrected leaves if the
    `corrections_summary` dict is provided. The output tree is written at
    outfolder/prefix_treeid.nhx.

    Args:
        tree (ete3.Tree): input tree to reconcile.
        outfolder (str): path to write the output
        tree_id (str): identifier of the tree, used in the output .nhx file name.
        sptree (str): name of the species tree file
        ali (str, optional): the fasta multiple alignment, required if branch lengths have to be
                             computed
        prefix (str, optional): string to add as prefix to the output file
        corrections_summary (dict, optional): dict of list of 5 element tuples summarizing
                                              corrections, the last tuple element should be a dict
                                              containing features to add to each leaf as .nhx tags.
        brlengths (bool, optional): Whether branch-lengths should be computed

    Returns:
        bool: True if no Exception was raised.
    """

    try:

        ete3_format = 1
        if not brlengths:
            ete3_format = 9

        out_exist = os.path.exists(outfolder+"/"+prefix+"_"+treeid) and\
                    os.path.getsize(outfolder+"/"+prefix+"_"+treeid) > 0

        if not resume or not out_exist:

            if brlengths:
                sys.stderr.write("Computing branch lengths for tree number "+treeid+"\n")
                sys.stderr.flush()

            wtree = Tree(tree, format=1)
            d_sp = {}

            leaves = wtree.get_leaves()

            #remove artefactual single-child nodes
            wtree.prune(leaves)

            #add species tag for treebest
            for leaf in leaves:
                d_sp[leaf.name] = leaf.S
                leaf.name = leaf.name +'_'+leaf.S

            # if requested, we re-compute branch lengths
            if brlengths:

                #extract ali and tree with species tag
                seq = ut.get_subali(ali, d_sp.keys(), d_sp)
                ut.write_fasta(seq, outfolder+"/tmp_"+treeid+".fa")
                wtree.write(outfile=outfolder+"/tmp_"+treeid, format=1, features=["S"],
                            format_root_node=True)

                #compute branch-length
                os.system("treebest phyml -t opt -n "+outfolder+"/tmp_"+treeid+".fa "+\
                          outfolder+"/tmp_"+treeid+" -c 2 > "+outfolder+"/"+treeid)

                #remove temp
                os.remove(outfolder+"/tmp_"+treeid+".fa")
                os.remove(outfolder+"/tmp_"+treeid)

            #otherwise, we just write the tree to file
            else:

                wtree.write(outfile=outfolder+"/"+treeid, format=ete3_format,
                            format_root_node=True)

            #Reconcile the tree
            os.system("treebest sdi -s "+sptree+" "+outfolder+"/"+treeid+\
                      " > "+outfolder+"/"+prefix+"_"+treeid)

            #remove temp
            os.remove(outfolder+"/"+treeid)

        wtree = Tree(outfolder+"/"+prefix+"_"+treeid, format=1)
        for leaf in wtree.get_leaves():
            leaf.name = leaf.name.replace('_'+leaf.S, '', 1)

        edit_tags = []

        if corrections:

            #extract dictionary storing a summary of modification for each leaf
            all_corrections_features = [d_feat for _, _, _, _, d_feat in corrections[int(treeid)]]

            #sorry about this ugly one-liner:
            #store for each leaf, WGD(s) for which it was corrected or had to be moved in the tree
            all_corrections_features = {key:list(itertools.chain(*[d_feat[key]\
                                        for d_feat in all_corrections_features\
                                        if key in d_feat])) for key in {key for d_feat in\
                                        all_corrections_features for key in d_feat}}
            edit_tags = gt.add_nhx_tags(wtree, all_corrections_features)

        #write tree
        all_features = ["S", "D", "DD", "DCS"] + edit_tags
        wtree.write(outfile=outfolder+"/"+prefix+"_"+treeid, format=ete3_format,
                    features=all_features, format_root_node=True)

        return True

    except Exception:

        traceback.print_exc()
        raise


def multiprocess_rec_brlgth(trees, alis, ncores, modified_trees, folder_cor, sptree,
                            prefix="cor", brlengths=True, resume=False):
    """
    Reconciles with the species tree and optionaly compute branch-lengths for a subset of trees in
    `modified_trees` of a gene trees forest, in parallel. Each output reconciled tree is written
    at outfolder/prefix_treeid.nhx.

    Args:
        trees (str): path to gene trees forest
        alis (str): path to the corresponding multiple alignments
        ncores (int): number of cores to use for parallel execution
        modified_trees (dict): gene trees to reconcile
        folder_cor (str): path to store each output reconciled tree file
        sptree (str): name of the species tree file
        prefix (str, optional): prefix to add to output trees
        brlengths (bool, optional): Whether branch-lengths should be computed

    """

    pool = multiprocessing.Pool(ncores)

    open_f = open
    if alis.split('.')[-1] == 'gz':
        open_f = gzip.open

    with open(trees, "r") as infile_t, open_f(alis, "rt") as infile_a:

        async_res = []

        for i, (input_tree, input_ali) in enumerate(zip(ut.read_multiple_objects(infile_t),
                                                        ut.read_multiple_objects(infile_a))):

            if i in modified_trees:
                res = pool.apply_async(worker_rec_brlgth, args=(input_tree, folder_cor, str(i),
                                                                sptree, input_ali, prefix,
                                                                modified_trees, brlengths, resume))
                async_res += [res]

        pool.close()
        pool.join()

    for res in async_res:
        if not res.get():
            sys.stderr.write("An error occured in a child process\n")
            sys.exit(1)



if __name__ == '__main__':

    #positional

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--treesFile', help='Forest of trees in the nhx format with species,\
                         duplication/speciation nodes + duplication confidence tags.',
                        required=True)

    PARSER.add_argument('-a', '--alisFile', help='Multiple fasta containing all multiple alignment\
                        used to build the gene trees of the forest.', required=True)

    PARSER.add_argument('-s', '--Species_tree', help='Newick species tree with ancestors',
                        required=True)

    PARSER.add_argument('-acc', '--Accepted',
                        help='File with subtrees to correct, accepted by the AU-test.',
                        required=True)

    PARSER.add_argument('-o', '--out', help='outfile name for the corrected forest',
                        required=True)


    PARSER.add_argument('-ogr', '--outgroup', help='outgroup_species', required=True)

    PARSER.add_argument('-anc', '--anc_dup', help='anc of duplicated_sp', required=True)

    #optional

    PARSER.add_argument('-tmp', '--tmpfolder',
                        help='folder for corrected trees with correction tags for visualization.',
                        required=False, default='tmp')

    PARSER.add_argument('-n', '--ncores', help='number of cores to use', required=False, default=1)


    PARSER.add_argument('-sa', '--save_cor', help='outfile name for the corrected forest',
                        required=False, default='n')

    PARSER.add_argument('-br', '--branch_lengths', help='Should branch lengths be computed ?',
                        required=False, default='y')

    PARSER.add_argument('-res', '--resume', help='Resume an interrupted run ?',
                        required=False, default='n')

    ARGS = vars(PARSER.parse_args())

    for ARG in ['save_cor', 'branch_lengths', 'resume']:
        assert ARGS[ARG] in ['y', 'n'], '{} should be y or n'.format(ARG)

        if ARGS[ARG] == 'y':
            ARGS[ARG] = True
        else:
            ARGS[ARG] = False


    #if it does not exist, create folder to store all corrected trees for visualization
    CORFOLDER = ARGS["tmpfolder"]

    if CORFOLDER:

        try:
            os.mkdir(CORFOLDER)

        except OSError:
            pass


    #Number of cores to use
    NCORES = int(ARGS["ncores"])

    #All ancestors with a WGD in the dataset
    WGD_ANCS = ARGS["anc_dup"].split(',')

    #Corresponding outgroup species used in the synteny analysis
    OUTGROUPS = ARGS["outgroup"].split('_')

    assert len(OUTGROUPS) == len(WGD_ANCS), "inconsistent ARGS for multiple wgds"

    #Load the list of corrected subtrees
    CORRECTED_SUBTREES = gt.load_corrections(ARGS['Accepted'])

    CORRECTION_STATS = {}
    WGDS = spt.get_anc_order(ARGS["Species_tree"], WGD_ANCS, tips_to_root=True)

    #Correct trees for each WGDs iteratively, from tips to root
    for j, wgd in enumerate(WGDS):


        TO_PRINT = (f"Re-grafting corrected subtrees for WGD {wgd}, outgroup species:"
                    f"{OUTGROUPS[WGD_ANCS.index(wgd)]}\n")
        sys.stderr.write(TO_PRINT)

        #start from last WGD-corrected forest (or input if j==0 i.e. first corrected WGD)
        if j == 0:
            TREES = ARGS["treesFile"]
        else:
            TREES = ARGS["out"]+'_'+str(j-1)

        #subtree corrected for this WGD
        CORRECTED_SUBTREES_CURRENT = {i:CORRECTED_SUBTREES[i] for i in CORRECTED_SUBTREES\
                                      if CORRECTED_SUBTREES[i][1] == wgd}

        #get list of species that are sister species of outgroups
        SISTERS_OF_OUTGR = {}
        for outg in OUTGROUPS[WGD_ANCS.index(wgd)].split(","):
            SISTERS_OF_OUTGR[outg] = spt.get_sister_species(ARGS['Species_tree'],\
                                                                   outg, wgd)
        #species below the current WGD only
        sp_wgd = spt.get_species(ARGS["Species_tree"], wgd, ARGS["anc_dup"])

        #get a list of species that underwent a subsequent wgd, as well as all their outgroups
        sp_other_wgd = {}
        for SUBS_WGD in WGDS[wgd]:
            sp_other_wgd[SUBS_WGD] = spt.get_species(ARGS["Species_tree"], SUBS_WGD)
            for sp in sp_wgd:
                SISTERS_OF_OUTGR[SUBS_WGD] = SISTERS_OF_OUTGR.get(SUBS_WGD, {})
                SISTERS_OF_OUTGR[SUBS_WGD][sp] = spt.get_sister_species(ARGS['Species_tree'],
                                                                        sp, SUBS_WGD)

        #correct the forest for the current WGD (topology only)
        with open(TREES, "r") as infile:

            for TREE_IND, TREE in enumerate(ut.read_multiple_objects(infile)):

                correct_wtrees(TREE, CORRECTED_SUBTREES_CURRENT, CORRECTION_STATS,
                               TREE_IND, CORFOLDER+'/cor_', SISTERS_OF_OUTGR,
                               sp_other_wgd, sp_wgd, tag=wgd)

        #Write the new forest and print some statistics
        ut.write_forest(TREES, ARGS["out"]+'_'+str(j), CORRECTION_STATS, wgd)

        #clean temp
        if j != 0:
            os.remove(TREES)


    j = len(WGDS) - 1

    #Now topological correction for all WGDs have been made
    sys.stderr.write('All topological corrections to the gene trees have been applied\n')

    if ARGS["branch_lengths"]:
        TO_PRINT = ("Reconciling corrected trees with the species tree & recomputing branch-"
                    "lengths\n")
        sys.stderr.write(TO_PRINT)
    else:
        sys.stderr.write('Reconciling corrected trees with the species tree\n')

    #reconcile with the species tree and re-compute branch-lengths
    multiprocess_rec_brlgth(ARGS["out"]+'_'+str(j), ARGS["alisFile"], NCORES, CORRECTION_STATS,
                            CORFOLDER, ARGS["Species_tree"], brlengths=ARGS["branch_lengths"],
                            resume=ARGS["resume"])

    #clean temp
    os.remove(ARGS["out"]+'_'+str(j))

    sys.stderr.write('Writing SCORPiOs-corrected gene trees forest\n')

    #write final forest
    ut.write_forest(ARGS["treesFile"], ARGS["out"], CORRECTION_STATS,
                    save_single_treefile=ARGS["save_cor"], cor_treefiles=CORFOLDER+'/cor_')

    #safely delete empty tmp_corrections folder if individual corrections are not to be saved
    if not ARGS["save_cor"]:

        #check that directory is empty
        if not os.listdir(CORFOLDER):

            os.rmdir(CORFOLDER)
