"""
This scripts allows to recover correction tags for trees that have been corrected multiple times
during iterative correction (currently .nhx tags are wiped out if a same tree is corrected again)

Example:
    $ python -m scripts.trees.iteration_nhx_tags -i 5
                                                 -c SCORPiOs_example/corrected_forest_%d.nhx
                                                 [-o out.nhx]
"""

import sys
import argparse

from ete3 import Tree

from . import utilities as ut


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--iter', help='Total number of iterations', required=True, type=int)

    PARSER.add_argument('-c', '--cor_f', help='path to corrected forests', required=True, type=str)

    PARSER.add_argument('-o', '--out', type=str, required=False, default="out.nhx")

    ARGS = vars(PARSER.parse_args())

    CORRECTIONS = {}

    COR_FOREST = ARGS["cor_f"]

    COR_TAGS_ALL = []

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

                for i in t.get_leaves():

                    corr = [att for att in vars(i) if "CORR_ID_" in att]

                    for tag in corr:

                        CORRECTIONS[i.name] = CORRECTIONS.get(i.name, [])

                        CORRECTIONS[i.name].append((tag+'_'+str(itera), getattr(i, tag)))

                        if tag+'_'+str(itera) not in COR_TAGS_ALL:

                            COR_TAGS_ALL.append(tag+'_'+str(itera))

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

                    tree_nhx = t.write(format=1, features=all_features, format_root_node=True)

                    OUTFILE.write(tree_nhx+'\n//\n')
