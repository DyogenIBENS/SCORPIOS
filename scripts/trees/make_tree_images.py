
"""
This script generates figures of SCORPiOs corrected trees and their uncorrected counterpart. It is
specifically designed to visualize SCORPiOs subtree corrections. Therefore, it assumes that inputs
are SCORPiOs-generated files, with SCORPiOs file naming and format conventions.

Leaves of each SCORPiOs-corrected wgd subtree are printed with the same color in both original and
corrected tree. Optionally, `--show_moved` also assigns a matching lighter color to leaves of
non-wgd species that have been rearranged to reinsert the wgd subtree.

Input can either be a list of files, in any order, containing any number of corrected/uncorrected
tree pairs, or a directory. The name of the corrected wgd and of the outgroups used should also be
provided.

Examples:

    $ python scripts/trees/make_tree_images.py \
-i SCORPiOs_example/Corrections/tmp_whole_trees_0/cor_27 \
SCORPiOs_example/Corrections/tmp_whole_trees_0/ori_27 --wgd Salmonidae --outgr \
Gasterosteus.aculeatus

    $ python scripts/trees/make_tree_images.py -i SCORPiOs_example/Corrections/tmp_whole_trees_0/ \
--wgd Clupeocephala --outgr Lepisosteus.oculatus,Amia.calva -f pdf -o img_clup --color_outgr

"""

import sys
import os
import argparse
import glob

from ete3 import Tree, NodeStyle, TreeStyle, TextFace


def color_internal_node(node, is_corrected_wgd=False):

    """
    Colors an internal node with convention colors: red for duplication, blue for speciation, cyan
    for dubious duplication.

    Arg:
        node (ete3.TreeNode): node to color
        is_corrected_wgd (bool, optional): set special style if node is corrected wgd node
    """

    nstyle = NodeStyle()

    nstyle["size"] = 4

    if hasattr(node, "D") and hasattr(node, "DD"):
        if node.DD == 'Y':
            nstyle["fgcolor"] = "cyan"

        elif node.D == 'Y':
            nstyle["fgcolor"] = "red"

    elif hasattr(node, "D"):
        if node.D == 'Y':
            nstyle["fgcolor"] = "red"

    elif hasattr(node, "DD"):
        if node.DD == 'Y':
            nstyle["fgcolor"] = "cyan"

    else:
        nstyle["fgcolor"] = "blue"

    #Annotate corrected wgd node with a big square.
    if is_corrected_wgd:
        nstyle["size"] = 10
        nstyle["shape"] = "sphere"
        nstyle["bgcolor"] = "#F0F0F0"

    node.set_style(nstyle)


def identify_outgroup_vs_wgd_subtree(subtrees, outgroups, tree):

    """
    Identifies outgroup subtree amongst given trees.

    Args:
        subtrees (list of ete3.TreeNode): input subtrees
        outgroups (list of str): name of outgroups
        tree (ete3.Tree): whole tree

    Returns:
        ete3.TreeNode: the node corresponding to wgd corrected subtree
    """
    for subtree in subtrees:
        leaves = subtree.get_leaves()
        if len(leaves) == 1:
            leaf = leaves[0]
            if leaf.S in outgroups:
                subtree_outgr = subtree
                break

    corrected_subtrees = {i for i in subtrees if i != subtree_outgr}
    if len(corrected_subtrees) == 1:
        corrected_subtree, = corrected_subtrees
    else:
        corrected_subtree = tree.get_common_ancestor(corrected_subtrees)
    return corrected_subtree


def get_corrected_wgd_nodes(tree, wgd, outgroups):

    """
    Finds all nodes that correspond to corrected wgd nodes.

    Args:
        tree (ete3.Tree): input tree
        wgd (str): restricts search to specific wgd

    Returns:
        list: the list of matched ete3.TreeNodes
    """

    nodes = []
    leaves = {i for i in tree.get_leaves()}
    all_attr = {getattr(i, "CORR_ID_"+wgd) for i in leaves if hasattr(i, "CORR_ID_"+wgd)}

    for attr in all_attr:
        corrected_subtrees = list(tree.get_monophyletic(values=[attr], target_attr="CORR_ID_"+wgd))
        assert corrected_subtrees
        if len(corrected_subtrees) == 1:
            corrected_subtrees = corrected_subtrees[0].get_children()

        corrected_wgd_node = identify_outgroup_vs_wgd_subtree(corrected_subtrees, outgroups, tree)
        nodes.append(corrected_wgd_node)

    return nodes


def color_leaves(node, palette, edit_d, wgd, usedict=False, moved=False, ignore=False):

    """
    Colors names of leaves of corrected subtrees.

    Args:
        node (ete3.TreeNode): leaf to color
        palette (list): pre-defined list of colors to use
        edit_d (dict): dictionary to store/load colors
        wgd (str): restrict coloring to specified wgd
        usedict (bool, optional): should dictionary be used to load colors
        moved (bool, optional): color rearranged non-wgd species leaves
    """

    #Hide terminal node
    node.img_style["size"] = 0
    color = 'black'

    #use nhx correction tags
    if not usedict:

        if not ignore:
            if hasattr(node, "CORR_ID_"+wgd):
                corr_id = (int(getattr(node, "CORR_ID_"+wgd)) - 1) % (len(palette)/2)
                color = palette[int(corr_id)*2]
                edit_d[node.name] = color

            elif moved and hasattr(node, "MOVED_ID_"+wgd):
                moved_id = (int(getattr(node, "MOVED_ID_"+wgd)) - 1) % (len(palette)/2)
                color = palette[int(moved_id)*2+1]
                edit_d[node.name] = color

    #use stored colors
    else:
        if node.name in edit_d:
            color = edit_d[node.name]

    name_face = TextFace(node.name, fgcolor=color, fsize=10)

    node.add_face(name_face, column=0, position='branch-right')


def make_tree_figure(tree, palette, outfile, edit_d, wgd, outgroups, usedict=False, outformat="png",
                     moved=False, coutgr=False):

    """
    Creates and saves an image for a tree object.

    Args:
        tree (ete3.Tree): input tree
        palette (list): pre-defined list of colors to use
        outfile (str): name for the output image
        edit_d (dict): dictionary to store/load colors
        wgd (str): restrict coloring to specified wgd
        usedict (bool, optional): should dictionary be used to load colors
        outformat (str, optional): output format for image png (default), svg or pdf
        moved (bool, optional): color rearranged non-wgd species leaves

    Returns:
        dict: colors to leaf names mapping
    """

    #create a treestyle that we'll customize
    tstyle = TreeStyle()
    tstyle.show_leaf_name = True

    corrected_wgd_nodes = get_corrected_wgd_nodes(tree, wgd, outgroups)

    if not corrected_wgd_nodes and not usedict:
        sys.stderr.write(f"No corrected subtree for wgd {wgd}, no image created.\n")
        return edit_d

    leaves = {i.name for i in tree.get_leaves()}
    if usedict and not leaves.intersection(set(edit_d.keys())):
        sys.stderr.write(f"No corrected subtree for wgd {wgd}, no image created.\n")
        return edit_d

    for node in tree.traverse():

        if not node.is_leaf():

            is_corrected_wgd = node in corrected_wgd_nodes

            #color duplication, speciation and dubious nodes with convention colors
            color_internal_node(node, is_corrected_wgd)

        #color edited leaves
        else:
            ignore = False
            if not coutgr and node.S in outgroups:
                ignore = True
            color_leaves(node, palette, edit_d, wgd, usedict=usedict, moved=moved, ignore=ignore)

    #color root node
    color_internal_node(tree)

    #hide default leaf name since we added custom ones
    tstyle.show_leaf_name = False

    #set width
    tstyle.tree_width = 200

    #save image
    out = outfile+"."+outformat
    tree.render(out, dpi=80, tree_style=tstyle)
    sys.stderr.write(f"{out} created !!\n")

    return edit_d


def make_all_figures(files, palette, edit_d, outfolder, wgd, outgrs, usedict=False, outformat="png",
                     moved=False, coutgr=False):

    """
    Generates images for a list of input files.

    Args:
        files (list): list of input files
        palette (list): pre-defined list of colors to use
        edit_d (dict): dictionary to store/load colors
        wgd (str): restrict coloring to specified wgd
        usedict (bool, optional): should dictionary be used to load colors
        outformat (str, optional): output format for image png (default), svg or pdf
        moved (bool, optional): color rearranged non-wgd species leaves


    Returns:
        dict: colors to leaf names mapping
    """

    for treef in files:
        sys.stderr.write(f"Reading {treef}...\n")

        outfile = outfolder + '/img_'+ os.path.splitext(os.path.basename(treef))[0]
        tree = Tree(treef)
        edit_d = make_tree_figure(tree, palette, outfile, edit_d, wgd, outgrs, usedict=usedict,
                                  outformat=outformat, moved=moved, coutgr=coutgr)
        sys.stderr.write("\n")

    return edit_d


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    REQUIRED = PARSER.add_argument_group('required named arguments')

    REQUIRED.add_argument('-i', '--input', nargs='+', help="Folder with corrected and original "
                          "trees, or a list of tree files.", required=True)

    REQUIRED.add_argument('-w', '--wgd', help='Corrected wgd to highlight. For instance, '
                          '-wgd Clupeocephala will show only subtrees corrected for the wgd that '
                          'occured in the Clupeocephala ancestor.', required=True)

    REQUIRED.add_argument('--outgr', help='Outgroup(s) used in SCORPiOs tree correction, '
                          'comma-separated.', required=True)


    PARSER.add_argument('-o', '--output', help='Output folder, default is trees_img/',
                        required=False, default='trees_img')


    PARSER.add_argument('-f', '--format', help='Output format (pdf, svg or png).', required=False,
                        default='png')

    PARSER.add_argument('--show_moved', help='Color non-wgd rearranged leaves, default is False',
                        action='store_true')

    PARSER.add_argument('--color_outgr', help='Color the outgroup gene used by SCORPiOs,'
                        'default is False', action='store_true')

    ARGS = vars(PARSER.parse_args())


    ##Check args

    #Check format
    ARGS["format"] = ARGS["format"].strip('.')
    assert ARGS["format"] in ['png', 'svg', 'pdf'], "output format should be pdf, svg or png\n"


    #Check input if it is a directory
    if len(ARGS["input"]) == 1 and not os.path.isdir(ARGS["input"][0]):
        INPUT = ARGS["input"][0]
        ERR = f"The provided input {INPUT} is not an existing directory, please check.\n"
        sys.stderr.write(ERR)
        PARSER.print_help(sys.stderr)
        sys.exit(1)


    #Load input if is a directory
    if len(ARGS["input"]) == 1:
        INFOLDER = ARGS["input"][0]
        COR = glob.glob(INFOLDER+'/cor_*')
        ORI = glob.glob(INFOLDER+'/ori_*')


    #Check input argument if it is list of files
    if len(ARGS["input"]) > 1 and len(ARGS["input"])%2 != 0:
        ERR = "You provided an odd number of input, please check.\n"
        sys.stderr.write(ERR)
        PARSER.print_help(sys.stderr)
        sys.exit(1)


    #Check and load input if it is list of files
    if len(ARGS["input"]) > 1:
        ORI, COR = [], []
        for input_file in ARGS["input"]:
            if os.path.basename(input_file)[0:3] == 'cor':
                COR.append(input_file)
                original_tree = input_file.replace("/cor_", "/ori_", 1)
                if original_tree not in ARGS["input"]:
                    sys.stderr.write(f"Missing original tree for tree {input_file}\n")
                    PARSER.print_help(sys.stderr)
                    sys.exit(1)

            elif os.path.basename(input_file)[0:3] == 'ori':
                ORI.append(input_file)
                cor_tree = input_file.replace("/ori_", "/cor_", 1)
                if cor_tree not in ARGS["input"]:
                    sys.stderr.write(f"Missing corrected tree for tree {input_file}\n")
                    PARSER.print_help(sys.stderr)
                    sys.exit(1)


    #define color palette
    PALETTE = ["crimson", "lightcoral", "steelblue", "lightblue", "darkgreen", "mediumseagreen",
               "brown", "peru", "mediumorchid", "plum", "deeppink", "hotpink", "lime",
               "greenyellow", "orangered", "darkorange", "goldenrod", "yellow"]


    #create output directory
    OUTFOLDER = ARGS["output"]
    os.makedirs(OUTFOLDER, exist_ok=True)

    #outgroups
    OUTGR = ARGS["outgr"].split(',')

    ## Create images

    COLORS_DICT = {}

    #create images for corrected tree and store used colors in a dictionary
    COLORS_DICT = make_all_figures(COR, PALETTE, COLORS_DICT, OUTFOLDER, ARGS["wgd"], OUTGR,
                                   usedict=False, outformat=ARGS["format"],
                                   moved=ARGS["show_moved"], coutgr=ARGS["color_outgr"])

    #create images for original tree using stored colors
    _ = make_all_figures(ORI, PALETTE, COLORS_DICT, OUTFOLDER, ARGS["wgd"], OUTGR, usedict=True,
                         outformat=ARGS["format"], moved=ARGS["show_moved"])
