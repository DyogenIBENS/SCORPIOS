"""
This scripts allows to recover correction tags for trees that have been corrected multiple times
during iterative correction (currently .nhx tags are wiped out if a same tree is corrected again)
Optionally also adds tag to internal corrected nodes that are corrected subtrees.

Example::

    $ python -m scripts.trees.iteration_nhx_tags -i 5
    -c SCORPiOs_example/corrected_forest_%d.nhx [-o out.nhx] [--internal]
"""

import sys
import argparse

from ete3 import Tree

from scripts.trees import speciestree as spt, utilities as ut
from scripts.synteny.duplicated_families import tag_duplicated_species


def corr_tag_below_node(node, tags_corr):

    """
    Search for the presence of .nhx tags for leaves below the node `node`.

    Args:
        node (ete3 TreeNode): the input node
        tags_corr (list of str): list of tags to search for

    Returns:
        bool: True if at least one of the input `tags_corr` is in leaves below `node`
    """

    has_tag = False

    for leaf in node.get_leaves():
        for tag_corr in tags_corr:
            if hasattr(leaf, tag_corr):
                has_tag = True
                break

    return has_tag


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--iter', help='Total number of iterations', required=True, type=int)

    PARSER.add_argument('-c', '--cor_f', help='path to corrected forests', required=True, type=str)

    PARSER.add_argument('--internal', help='tag also internal corrected wgd nodes',
                        action='store_true')

    PARSER.add_argument('-sp', '--sptree', type=str, required=False, default='')

    PARSER.add_argument('-o', '--out', type=str, required=False, default="out.nhx")

    ARGS = vars(PARSER.parse_args())

    assert (ARGS["internal"] and ARGS["sptree"]) or not ARGS["internal"],\
           "A species tree should be provided to tag corrected wgd internal nodes"

    CORRECTIONS = {}

    COR_FOREST = ARGS["cor_f"]

    COR_TAGS_ALL = []

    COR_TAGS_INT = []

    COR_TAGS = []


    for itera in range(1, ARGS["iter"] + 1):
        k = 0

        if itera != ARGS["iter"]:

            sys.stderr.write(f"Browsing corrected forest at iteration {itera}\n")

        else:
            eprint = (f"Browsing corrected forest at final iteration ({itera}) and writing final "
                      f"forest with tags\n")
            sys.stderr.write(eprint)

        with open(COR_FOREST % itera, 'r') as f, open(ARGS["out"], 'w') as OUTFILE:

            for tree in ut.read_multiple_objects(f):

                if k%1000 == 0 and k > 0:
                    sys.stderr.write(f"Browsed {k} trees\n")

                k += 1

                t = Tree(tree)
                leaves = t.get_leaves()

                for i in leaves:

                    corr = [att for att in vars(i) if "CORR_ID_" in att]

                    for tag in corr:

                        CORRECTIONS[i.name] = CORRECTIONS.get(i.name, [])

                        CORRECTIONS[i.name].append((tag+'_'+str(itera), getattr(i, tag)))

                        if tag+'_'+str(itera) not in COR_TAGS_ALL:

                            COR_TAGS_ALL.append(tag+'_'+str(itera))

                        if tag not in COR_TAGS_INT:
                            COR_TAGS_INT.append(tag)

                    if itera == ARGS["iter"]:
                        if i.name in CORRECTIONS:
                            seen = []
                            for tag, value in CORRECTIONS[i.name]:
                                wgd = tag.split("_")[-2]
                                if wgd not in seen:
                                    setattr(i, tag, value)
                                    seen.append(wgd)

                if itera == ARGS["iter"]:

                    all_features = ["S", "D", "DD", "DCS"] + COR_TAGS_ALL

                    if not ARGS["internal"]:

                        tree_nhx = t.write(format=1, features=all_features, format_root_node=True)

                        OUTFILE.write(tree_nhx+'\n//\n')


                    else:

                        wgds = {i.split('_')[-2] for i in COR_TAGS_ALL}
                        wgds_dict = {}
                        for wgd in wgds:
                            DUP_SPECIES = spt.get_species(ARGS["sptree"], wgd)
                            wgds_dict[wgd] = DUP_SPECIES

                        t = Tree(t.write(format=1, features=all_features, format_root_node=True))
                        for wgd in wgds_dict:

                            leaves = t.get_leaves()
                            if len(leaves) == 1:
                                continue

                            #find all monphyletic teleost groups
                            tag_duplicated_species(leaves, wgds_dict[wgd])

                            #all clades of teleost genes,
                            #by definition corrected subtrees will only contain dup. sp
                            subtrees = t.get_monophyletic(values=["Y"], target_attr="duplicated")

                            for subtree in subtrees:

                                if subtree.is_leaf():
                                    continue

                                #if corrected leaves at each side of the node: corrected node
                                child1, child2 = subtree.get_children()

                                tags_wgd = [i for i in COR_TAGS_ALL if wgd in i]

                                ok_child1 = corr_tag_below_node(child1, tags_wgd)

                                if ok_child1:
                                    ok_child2 = corr_tag_below_node(child2, tags_wgd)

                                if ok_child1 and ok_child2:

                                    setattr(subtree, 'CORR_ID_'+wgd, 'Y')

                            t = Tree(t.write(format=1, features=all_features+COR_TAGS_INT,
                                             format_root_node=True))

                        tree_nhx = t.write(format=1, features=all_features+COR_TAGS_INT,
                                           format_root_node=True)

                        OUTFILE.write(tree_nhx+'\n//\n')
