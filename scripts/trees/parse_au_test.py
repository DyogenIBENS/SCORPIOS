#!/usr/bin/env python


"""
    Script to parse results of gene trees likelihood AU-test (here for comparison of 2 trees only),
    as written by CONSEL (Shimodaira, 2002).

    Example::

        $ python -m scripts.trees.parse_au_test -i inputs_polyS.txt [-o Accepted_Trees]
        [-it inputs_treeB.txt] [-one n] [-wgd Clupeocephala] [-p path/tree.nh] [--lore]
"""

import os
import sys
import argparse

#FIXME: DRY once lorelei mode completly implemented... (low priority):
# Write only one function to parse CONSEL results, regardless of the number of trees

def one_file_consel(filename, alpha, item_test='1'):

    """
    Parses one consel file record to obtain AU-test results. This function assumes that the
    AU-test was performed to compare two trees only.

    By default, the tree labelled '1' in consel is considered the tested tree, result indicate
    whether this tree is as good as (or better) as the other tree (reference).

    Args:
        fi (str): input filename.
        alpha (float, optional): alpha threshold for significance of the AU-test.
        item (str, optional): tested tree consel label, the other is considered the reference.

    Returns:
        str: One of 'error', 'rejected', 'equivalent lower lk', 'equivalent higher lk' or
        'better sign. higher lk'.

    Note:

        The function returns the 'error' string in two cases:

        - the input file is empty, in SCORPiOs workflow this happens when the synteny aware tree
          could not be built with ProfileNJ due to fastdist failing to build the distance matrix.

        - CONSEL failed to compute the log-likelihoods from the phyml or raxml site likelihood file
          --> logs should be checked.
    """

    #if consel failed (probably because tree building failed) we return the 'error' value.
    au_result = 'error'

    if  os.stat(filename).st_size != 0:

        tmp = []

        with open(filename, 'r') as infile:

            for line in infile:

                if line[0] == '#' and line[0:3] != '# r':

                    res = line.strip().split()

                    #If no error occurred each result line in consel should contain more than 4 col
                    if len(res) > 4:

                        item = res[2]
                        au_value = res[4]
                        tmp.append((item, au_value))
                        obs = res[3]

            if tmp and obs != 'inf' and len(tmp) > 1:

                if tmp[0][0] == item_test:

                    au_result = 'equivalent higher lk'

                    if float(tmp[1][1]) < alpha:

                        au_result = 'better sign. higher lk'

                elif float(tmp[1][1]) >= alpha:

                    au_result = 'equivalent lower lk'

                else:

                    au_result = 'rejected'

            elif tmp:

                sys.stderr.write(f"Warning: for {filename}, CONSEL failed to compute log-lk\n")

    return au_result


def one_file_consel_3_trees(filename, alpha, item_dict=None):

    """
    Parses a CONSEL result file for a comparison of 3 trees.

    Args:
        filename (str) : name of CONSEL file to parse
        alpha (float): alpha threshold for significance of the AU-test.
        item_dict (dict, optional): correspondance between item in CONSEL and tree labels

    Returns:
        str: One of 'error', 'convergence_pb', 'aore rejected', 'lore rejected' or
        'lore and aore rejected'.
    """

    #if consel failed (probably because tree building failed) we return the 'error' value.
    au_result = 'error'

    #correspondance between item ids in CONSEL and tree types
    if item_dict is None:
        item_dict = {"ml": "1", "aore": "2", "lore":"3"}

    if  os.stat(filename).st_size != 0:

        tmp = []

        with open(filename, 'r') as infile:

            for line in infile:

                if line[0] == '#' and line[0:3] != '# r':

                    res = line.strip().split()

                    #If no error occurred each result line in consel should contain more than 4 col
                    if len(res) > 4:

                        item = res[2]
                        au_value = res[4]
                        tmp.append((item, au_value))
                        obs = res[3]

                        if obs != 'inf':
                            if item == item_dict['ml'] and float(au_value) < alpha:
                                au_result = 'convergence_pb'
                                break

                            elif item == item_dict['lore'] and float(au_value) < alpha:
                                if au_result == 'aore rejected':
                                    au_result = 'lore and aore rejected'
                                else:
                                    au_result = 'lore rejected'

                            elif item == item_dict['aore'] and float(au_value) < alpha:
                                if au_result == 'lore rejected':
                                    au_result = 'lore and aore rejected'
                                else:
                                    au_result = 'aore rejected'

            if tmp and obs != 'inf' and len(tmp) > 1:

                if au_result == 'error':

                    au_result = 'neither lore or aore rejected'

            elif tmp:

                sys.stderr.write(f"Warning: for {filename}, CONSEL failed to compute log-lk\n")

    return au_result


def count(filenames, name_sol="", alpha=0.05, item='1', parse_name=True, wgd=''):

    """
        Parses all consel outputs in the input list and returns a list of 'accepted trees'.
        Prints to screen the numbers and proportion of trees:
        (i) accepted by the likelihood test
        (ii) accepted with a better likelihood
        (iii) accepted with a significantly better likelihood.

        Args:
            filenames (list of str): List of consel result files.
            name_sol (str, optional): Tag for tested trees that will be printed with the results.
                                      For instance, it can be the program used to build tested
                                      trees.
            alpha (float): alpha threshold for significance of the AU-test.
            item (str, optional): tested tree consel label, the other is considered the reference.
            parse_name (bool, optional): parse filename (expects SCORPiOs naming pattern)
            wgd (str, optional): name of the corrected wgd, to print out with the result summary

        Returns:
            list: list of names of accepted trees.
    """

    a_list = []
    all_res = {}

    for filename in filenames:


        #naming parsing is specific to SCORPiOs inputs and outputs
        #parsing below allows to retrieve the outgroup gene used as name for gene families.
        if parse_name:
            name = os.path.splitext(os.path.basename(filename).split('Res_')[1])[0]
        else:
            name = filename

        au_result = one_file_consel(filename, alpha, item)

        all_res[au_result] = all_res.get(au_result, 0) + 1

        if au_result not in ['rejected', 'error']:
            a_list.append(name)


    tot = sum(all_res.values())
    err = all_res.get('error', 0)

    if tot - err != 0:

        higher_lk = all_res.get('better sign. higher lk', 0) +\
                    all_res.get('equivalent higher lk', 0)
        acc = higher_lk + all_res.get('equivalent lower lk', 0)
        better = all_res.get('better sign. higher lk', 0)

        acc_prop = str(round(acc/float(tot-err)*100, 1))
        higher_p = str(round(higher_lk/float(tot-err)*100, 1))
        better_prop = str(round(better/float(tot-err)*100, 1))

        reason = '??: Check logs.'

        if name_sol != "ProfileNJ":
            reason = "Likely reason: ProfileNJ solution already accepted" # should be checked

        print('\n')
        print("----------------------------------AU-TESTs----------------------------------")
        print(" Whole-genome duplication: {}".format(wgd))
        print("\n")
        print(" Likelihood-tests results for {} solutions.".format(name_sol))
        print(" Total : {} tested subtrees                                 ".format(tot-err))
        print(" Accepted correction (similar or better lk):  {} ({}%)         ".format(acc,\
                                                                                       acc_prop))
        print(" Accepted correction with higher lk:  {} ({}%)       ".format(higher_lk,\
                                                                             higher_p))
        print(" Accepted correction with sign. higher lk:  {} ({}%) ".format(better, better_prop))
        if err:
            print(" {} subtrees not tested ({})".format(err, reason))
        print("----------------------------------------------------------------------------")
        print("\n")

    else:
        print('\n')
        print("----------------------------------AU-TESTs----------------------------------")
        print(" Whole-genome duplication: {}".format(wgd))
        print("\n")
        print(" Likelihood-tests results for {} solutions.".format(name_sol))
        print(" Total : {} tested subtrees                                 ".format(tot-err))
        print("----------------------------------------------------------------------------")
        sys.stderr.write("None of consel result files could be parsed, check the inputs")

    return a_list


def lore_aore_summary(filenames, alpha=0.05, item_dict=None, parse_name=True, wgd=''):

    """
    Parses all consel outputs in the input list and returns a summary, telling, for each tree, if
    the lore or aore (or both) topologies can be rejected.
    Also prints statistics to stdout.

    Args:
        filenames (list of str): List of consel result files.
        alpha (float): alpha threshold for significance of the AU-test.
        item_dict (str, optional): tested tree consel label, the other is considered the reference.
        parse_name (bool, optional): parse filename (expects SCORPiOs naming pattern)
        wgd (str, optional): name of the corrected wgd, to print out with the result summary

    Returns:
        dict: dictionary with the AU-tests results summary.

    """

    res_dict = {}
    all_res = {}

    for filename in filenames:

        #naming parsing is specific to SCORPiOs inputs and outputs
        #parsing below allows to retrieve the outgroup gene used as name for gene families.
        if parse_name:
            name = os.path.splitext(os.path.basename(filename).split('Res_')[1])[0]
        else:
            name = filename

        au_result = one_file_consel_3_trees(filename, alpha, item_dict)
        res_dict[name] = au_result
        all_res[au_result] = all_res.get(au_result, 0) + 1

    tot = sum(all_res.values())
    err = all_res.get('error', 0)

    if tot - err != 0:

        both_rej = all_res.get('lore and aore rejected', 0)
        both_acc = all_res.get('neither lore or aore rejected', 0)
        lore = all_res.get('aore rejected', 0)
        aore = all_res.get('lore rejected', 0)
        issue = all_res.get('convergence_pb', 0)

        reason = '??: Check logs.'

        print('\n')
        print("----------------------------------AU-TESTs----------------------------------")
        print(" LORe vs AORe at speciation: {}".format(wgd))
        print("\n")
        print(" Total : {} tested subtrees                                 ".format(tot-err))
        print(" AORe topologies (LORe rejected):  {}        ".format(aore))
        print(" LORe topologies (AORe rejected):  {}       ".format(lore))
        print(" Ambiguous : {} ({} both rejected, {} none rejected, {} convergence issue)"\
              .format(both_rej+both_acc+issue, both_rej, both_acc, issue))
        if err:
            print(" {} subtrees not tested ({})".format(err, reason))
        print("----------------------------------------------------------------------------")
        print("\n")

    else:
        print('\n')
        print("----------------------------------AU-TESTs----------------------------------")
        print(" LORe vs AORe at speciation: {}".format(wgd))
        print("\n")
        print(" Total : {} tested subtrees                                 ".format(tot-err))
        print("----------------------------------------------------------------------------")
        sys.stderr.write("None of consel result files could be parsed, check the inputs")

    return res_dict

if __name__ == '__main__':


    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--input', type=str, help='Either a file listing all input\
                        filenames or a single input filename', required=True)

    PARSER.add_argument('-o', '--output', help='Output file', default='out', required=False)

    PARSER.add_argument('-one', '--one_file', help='Is the input a file listing inputs to parse\
                        or a single file to parse.', required=False, default='n')

    PARSER.add_argument('-it', '--input2', help='A second file listing input filenames to be\
                        parsed separetely, for instance trees build with another method.',
                        required=False)

    PARSER.add_argument('-w', '--wgd', help='', required=False, default='')

    PARSER.add_argument('-p', '--path', help='Path to corresponding tree',
                        required=False, default='')

    PARSER.add_argument('--lore', help='Call script in LORelEi mode: i.e to compare ml, lore and\
                        aore trees', required=False, action='store_true')

    ARGS = vars(PARSER.parse_args())

    assert ARGS["one_file"] in ['y', 'n'], 'one_file should be y or n'

    if ARGS["one_file"] == 'y':
        ARGS["one_file"] = True
    else:
        ARGS["one_file"] = False

    ACCEPTED_TREES_1, ACCEPTED_TREES_2 = [], []

    #if the script is called on one file
    if ARGS['one_file']:
        AU_RESULT = one_file_consel(ARGS['input'], 0.05, '1')
        if AU_RESULT not in ['rejected', 'error']:
            print("true")
        else:
            print("false")

    #if the script is called on all results
    elif not ARGS["lore"]:
        TREE_PATHS = ARGS["path"].split(',')
        with open(ARGS['input'], 'r') as INFILE:
            INPUTS = [LINE.strip() for LINE in INFILE]
        ACCEPTED_TREES_1 = count(INPUTS, name_sol="ProfileNJ", wgd=ARGS['wgd'])

        if ARGS['input2']:

            with open(ARGS['input2'], 'r') as INFILE:
                INPUTS = [LINE.strip() for LINE in INFILE]
            ACCEPTED_TREES_2 = count(INPUTS, name_sol="TreeBest", wgd=ARGS['wgd'])

        with open(ARGS['output'], "w") as OUTFILE:
            for i in ACCEPTED_TREES_1:
                OUTFILE.write(i+'\t'+TREE_PATHS[0]+'\t'+ARGS['wgd']+'\n')

            for i in ACCEPTED_TREES_2:
                if i not in ACCEPTED_TREES_1:
                    OUTFILE.write(i+'\t'+TREE_PATHS[1]+'\t'+ARGS['wgd']+'\n')

    else:
        with open(ARGS['input'], 'r') as INFILE:
            INPUTS = [LINE.strip() for LINE in INFILE]
        SUMMARY_TREES = lore_aore_summary(INPUTS, wgd=ARGS['wgd'])

        with open(ARGS['output'], "w") as OUTFILE:
            for i in SUMMARY_TREES:
                OUTFILE.write(i+'\t'+SUMMARY_TREES[i]+'\n')
