#!/usr/bin/env python

"""
    Module with functions to combine SCORPiOs synteny predictions (orthogroups in graphs) across
    multiple outgroups.
"""


import sys
import itertools


def read_summary(input_file):

    """
    Loads a tab-delimited file summarizing the number of edges that were cut in each graphs to
    find the orthogroups.

    Args:
        input_file (str): input file name.

    Returns:
        dict: for each graph (key), the number of cuts (value).
    """

    summary = {}

    with open(input_file, 'r') as infile:

        for i, line in enumerate(infile):

            if i != 0:

                line = line.strip().split('\t')
                try:
                    fam, nb_cuts = line[0], int(line[-1])
                    summary[fam] = nb_cuts
                except ValueError:
                    pass
    return summary


def map_families_across_outgr(graphs):

    """
    Extracts corresponding graphs across outgroups.

    Args:
        graphs (OrderedDict of dict): For each outgroup (key level1), orthogroups in graphs of
                                      each families (key level 2),
                                      represented by a FamilyOrthologies object.

    Returns:
        combin (list of list): For all graphs, lists of all graphs_ids in all outgroups.
        mapped_ids (dict): For each graph_id, store its position in the combin list and the name of
                           the outgroup species.
    """

    combin = []
    mapped_ids = {}
    j = 0

    for i, outgr1 in enumerate(graphs):

        for family_id in graphs[outgr1]:

            genes = set(graphs[outgr1][family_id].orthogroup_a +\
                        graphs[outgr1][family_id].orthogroup_b)

            #for each outgroup species pair
            for outgr2 in list(graphs.keys())[i+1:]:

                for family_id2 in graphs[outgr2]:

                    genes2 = set(graphs[outgr2][family_id2].orthogroup_a +\
                                 graphs[outgr2][family_id2].orthogroup_b)

                    if genes.isdisjoint(genes2):
                        continue

                    #is family_id2 already in a family ?
                    #if yes add current family_id to that family
                    if family_id2 in mapped_ids:
                        combin[mapped_ids[family_id2][0]].append(family_id)
                        mapped_ids[family_id] = (mapped_ids[family_id2][0], outgr1)


                    #same for family_id
                    if family_id in mapped_ids:
                        combin[mapped_ids[family_id][0]].append(family_id2)
                        mapped_ids[family_id2] = (mapped_ids[family_id][0], outgr2)


                    #else add the group
                    else:
                        combin.append([family_id, family_id2])
                        mapped_ids[family_id2] = (j, outgr2)
                        mapped_ids[family_id] = (j, outgr1)
                        j += 1

    return combin, mapped_ids


def choose_best_graph(m_fam, mapped_ids, all_graphs, final_graphs, summ):

    """
    Chooses orthogroup prediction where the families were the less aggregated (one outgroup can
    over-aggregate families if one ortholog was lost) and community detection was easiest.

    Args:
        m_fam (list): a group of graphs matched across outgroups (same gene family)

        mapped_ids (dict): For each graph_id, store its position in m_fam and the name of
                           the outgroup species.

        all_graphs (OrderedDict of dict): For each outgroup (key level1), orthogroups in graphs of
                                          each families (key level 2),
                                          represented by a FamilyOrthologies object.

        final_graphs (dict): Dict to store results

        summ (dict): For each graph (key), the number of cuts (value).

    Returns:
        graphs (dict): for each outgroup, the gene id(s) of the combined family
        selected_graph (list): list of selected gene id(s) (family outgroup with best prediciton)
    """

    graphs = {}
    selected_graph = []
    for outgr in all_graphs:
        genes_in_outgr = {gene for gene in m_fam if mapped_ids[gene][1] == outgr}
        if genes_in_outgr:
            graphs[outgr] = genes_in_outgr

    #less aggregated families
    max_subfam = max(graphs.items(), key=lambda x: len(x[1]))
    all_max = [i[0] for i in graphs.items() if len(i[1]) == len(max_subfam[1])]

    #if 1 is less aggregated than all others, choose it
    if len(all_max) == 1:
        selected_graph = graphs[all_max[0]]
        for family_id in graphs[all_max[0]]:
            final_graphs[family_id] = all_graphs[all_max[0]][family_id]


    #otherwise choose the one with less edges cut in graph
    else:
        cuts = {}
        for outgr in graphs:
            cuts[outgr] = sum([summ[outgr][gene] for gene in graphs[outgr]])

        min_cut = min(cuts.items(), key=lambda x: x[1])
        all_min = [i[0] for i in cuts.items() if i[1] == min_cut[1]]
        if len(all_min) == 1:
            selected_graph = graphs[all_min[0]]
            for family_id in graphs[all_min[0]]:
                final_graphs[family_id] = all_graphs[all_min[0]][family_id]

        #otherwise choose first outgroup in prioritized list
        else:
            candidates = [i for i in all_graphs.keys() if i in all_min]
            selected_graph = graphs[candidates[0]]
            for family_id in graphs[candidates[0]]:
                final_graphs[family_id] = all_graphs[candidates[0]][family_id]
    return graphs, selected_graph


def combine_outgroups(all_graphs, summ_files, out='out'):

    """
    Combines orthogroups predictions across all outgroups using best graphs.

    Args:
        all_graphs (OrderedDict of dict): For each outgroup (key level1), orthogroups in graphs of
                                          each families (key level 2),
                                          represented by a FamilyOrthologies object.

        summ_files (str): Comma-separated file names summarizing graphs cuts for each outgroup,
                          outgroups should be in the same order as in all_graphs.

        outfile (str, optional): File to write a summary of outgroup graphs selected

    Returns:
        final_graphs (dict): For each family, orthogroup predictions of chosen graph represented
                             by a FamilyOrthologies object.
    """

    final_graphs = {}
    summ = {}

    for i, outgr in enumerate(all_graphs):
        summ[outgr] = read_summary(summ_files.split(',')[i])

    combin, mapped_ids = map_families_across_outgr(all_graphs)

    with open(out, 'w') as outfile:

        outfile.write('\t'.join(all_graphs.keys())+'\tbest_graph\n')

        #families only in one outgroup
        for i, outgr in enumerate(all_graphs):
            combin1d = list(itertools.chain.from_iterable(combin))
            only_in_one = [i for i in all_graphs[outgr] if i not in combin1d]

            for family in only_in_one:
                final_graphs[family] = all_graphs[outgr][family]

                #write summary to file
                genes_outgr = ['' if ogr != outgr else family for ogr in all_graphs]
                outfile.write("\t".join(genes_outgr)+'\t'+family+'\n')


        #for others choose the best graph
        for matched_fam in combin:
            graphs, families = choose_best_graph(matched_fam, mapped_ids, all_graphs,
                                                 final_graphs, summ)
            genes_outgr = []

            #write summary to file
            for outgr in all_graphs:
                genes = ''
                if outgr in graphs:
                    genes = (",").join(graphs[outgr])
                genes_outgr.append(genes)

            outfile.write("\t".join(genes_outgr)+'\t'+(",").join(families)+'\n')

    return final_graphs


if __name__ == '__main__':
    sys.exit()
