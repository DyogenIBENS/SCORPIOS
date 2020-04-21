#!/usr/bin/env python

"""
    Script to build orthology graphs and detect communities in the graphs.
    For each gene family, the script searches for two orthologous gene communities that split from
    a whole genome duplication.

    Example:
        $ python -m scripts.graphs.orthogroups -i orthology_file.gz [-o out] [-w n] [-n 1]
                                                                    [-s Summary] [-ignSg y]
                                                                    [-wgd '']  [--spectral]
                                                                    [--verbose]
"""

import sys
import argparse
import traceback
import multiprocessing
import operator

from collections import Counter
import gzip

import networkx as nx
from sklearn.cluster import SpectralClustering

def species_name(gene):

    """
    Parses gene name to extract species name. Expects species name after the last '_'.

    Args:
        gene (str): Gene name.

    Returns:
        str: Species name

    """

    species = gene.split('_')[-1]

    return species


def load_line(line, use_weights=False):

    """
    Loads a single orthology, which will be used to update the graphs.

    Args:
        line (str): a single line of the orthology file.

    Returns:
        edges (list of tuples): list of orthology pairs.
        weights (list): list of corresponding weights.
        fam_id (str): Unique id of the gene family
    """

    edges = []
    weights = None

    if use_weights:
        weights = []

    ortho1, ortho2, weight, fam_id = line.strip().split('\t')

    for gene in ortho2.split(','):
        gene1, gene2 = sorted([ortho1, gene])
        edges.append((gene1, gene2))
        if use_weights:
            weights.append(float(weight))

    return edges, weights, fam_id


def lazy_load_pairwise_file(file_object, use_weights=False):

    """
    Loads orthologies for a gene family and build the graph, one gene family at a time.
    The input should be a tab-delimited file, with the following columns:
    ortho_gene1, ortho_gene2, orthology confidence, gene family ID.

    Args:
        file_object (File Object): python file object for the input file
        use_weights (bool, optional): whether weights should be used in the graphs

    Yields:
        fam (networkx graph): Orthology graph of the gene family
        prev_id (str): Unique id of the gene family
    """

    prev_id = ''
    fam = nx.Graph()
    all_edges = []
    while True:
        line = file_object.readline()

        #if last line yield family and end
        if not line:
            yield fam, prev_id
            break

        fam_id = line.strip().split('\t')[-1]

        #if new family yield the previous and store the new one
        if prev_id and fam_id != prev_id:
            fam.add_edges_from(all_edges)
            all_edges = []
            yield fam, prev_id

            edges, weights, fam_id = load_line(line, use_weights)
            if use_weights:
                edges = edges[0] + (weights,)
            else:
                all_edges += edges
            # print(all_edges)
            fam = nx.Graph()

            # update_networkx(fam, edges, weights)
            prev_id = fam_id


        #if same family keep storing orthologies
        else:
            edges, weights, fam_id = load_line(line, use_weights)
            if use_weights:
                edges = edges[0] + (weights,)
                all_edges.append((edges, weights))
            else:
                all_edges += edges
            # print(all_edges)
            # update_networkx(fam, edges, weights)
            prev_id = fam_id


def most_central_edge(graph):

    """
    Extracts the most central edge of a graph, taking weights into account.
    If all weights are equal, the most central edge is the edge with the highest unweighted
    betweenness centrality.

    Args:
        graph (networkx.Graph): input graph

    Returns:
        networkx.edge: The most central edge.

    """

    #compute betweenness centrality of all edges, with weights.
    centrality = nx.edge_betweenness_centrality(graph, weight='weight')
    edge = max(centrality, key=centrality.get)

    return edge


def contracted_nodes(graph, u, v):

    """
    Modifies the graph by contracting `u` and `v`. Node contraction identifies the two nodes as a
    single node incident to any edge that was incident to the original two nodes.
    The right node `v` will be merged into the node `u`, so only `u` will appear in the
    returned graph.

    Args:
        graph (networkx.Graph): Orthology graph
        u, v (str, str) : name of nodes to contract, must be in `graph`.

    Note:
        Adapted from https://www.bountysource.com/issues/46183711-contracted_nodes-with-weights-\
                     giving-different-answers-according-to-order-of-inputs

    """

    new_edges = []
    nodes_u = graph[u]

    #browse nodes connected to v
    for w, feature_dict in graph[v].items():

        edge_weight = feature_dict.get('weight', 1)

        if w not in nodes_u:
            new_edges.append((u, w, edge_weight))

        else:
            #weight of common edges is set to max weight
            max_weight = max(edge_weight, graph[u][w].get('weight', 1))
            graph[u][w]['weight'] = max_weight

    #remove all edges connected to v
    v_data = graph.node[v]
    graph.remove_node(v)

    #add contracted edges
    if new_edges:
        graph.add_weighted_edges_from(new_edges)

    #keep contraction history
    graph.node[u]['contraction'] = graph.node[u].get('contraction', {})
    graph.node[u]['contraction'][v] = v_data


def collapse_nodes(graph):

    """
    Collapses tandem duplicates into a single node. Tandem duplicates are genes of a species
    having the same edges in the graph.

    Args:
        graph (networkx.Graph): Orthology graph
    """

    #identify nodes with same neighbors in graph (tandem duplicates)
    collapse_dict = {}
    for node in graph:
        neighbors = frozenset(graph.neighbors(node))
        species = species_name(node)

        if species not in collapse_dict:
            collapse_dict[species] = {}

        if neighbors not in collapse_dict[species]:
            collapse_dict[species][neighbors] = []

        collapse_dict[species][neighbors].append(node)

    #collapse
    for species in collapse_dict:
        for collapse in collapse_dict[species].keys():

            #if at least 2 nodes have exact same neighbor, collapse them
            if len(collapse_dict[species][collapse]) > 1:
                representative_node = collapse_dict[species][collapse][0]

                for node in collapse_dict[species][collapse][1:]:
                    contracted_nodes(graph, representative_node, node)


def are_species_sep(partitions):

    """
    Checks whether genes (or place-holders loss of a duplicate) of each species are split in the
    two communities.

    Args:
        partitions (tuple): Genes in each graph community in tuples of tuple

    Returns:
        bool: The return value, True if genes of the same species are in two communities, False
              otherwise.

    """

    #Sort by order of length, to extract only the two main partitions.
    partitions = list(partitions)
    partitions.sort(key=len, reverse=True)
    total_species = partitions[0].union(partitions[1])

    #Filter out species with only one representative.
    all_sp = [species_name(i) for i in total_species]
    count = Counter(all_sp)
    to_filter = [k for k in count if count[k] == 1]

    #Find species unique to partition1 or unique to partition2.
    seen_part1 = {species_name(i) for i in partitions[0] if species_name(i) not in to_filter}
    seen_part2 = {species_name(i) for i in partitions[1] if species_name(i) not in to_filter}
    non_sep = len(seen_part1 - seen_part2) + len(seen_part2 - seen_part1)

    return non_sep


def min_cut(graph, spectral=False):

    """
    Detects two orthologous communities in the graph.

    Args:
        graph (networkx.Graph): Orthology graph
        spectral (bool, optional): Use spectral clustering instead of default Girvan-Newman
                                   (faster)
    """

    scores_cut = []
    scores_uncut = []

    components = [i for i in nx.connected_components(graph)]
    all_sp = [species_name(i) for i in graph.nodes()]
    count = Counter(all_sp)


    #if only two genes, tree topology will be the same regardless of cuts
    if len([i for i in graph.nodes() if i != 'None_'+species_name(i)]) == 2:
        algo = 'Too few genes'
        cut_edges_best = {}

    #if more than one component, partitions are already defined
    elif len(components) > 1:
        cut_edges_best = {}
        algo = 'No'


    #if only one gene per species and only one component only one group
    elif max(count.items(), key=operator.itemgetter(1))[0] == 1:
        cut_edges_best = {}
        algo = 'No'

    # else find communities 
    else:

        #use spectral clustering only, if spectral clustering is preferred
        if spectral:
            algo = "spectral"
            adj_mat = nx.to_numpy_matrix(graph)
            sc = SpectralClustering(2, affinity='precomputed', n_init=100)
            sc.fit(adj_mat)
            components = ({n for i, n in enumerate(graph.nodes()) if sc.labels_[i]==0},\
                          {n for i, n in enumerate(graph.nodes()) if sc.labels_[i]==1})
            cut_edges_best = {(u, v) for u in components[0]\
                              for v in components[1].intersection(graph.adj[u])}

        #Default with girvan_newman and kerningan-lin
        else:
            comp = nx.algorithms.community.girvan_newman(graph, most_valuable_edge=most_central_edge)

            components = tuple(set(c) for c in next(comp))

            #edges between two sets of nodes
            cut_edges_best = {(u, v) for u in components[0]\
                              for v in components[1].intersection(graph.adj[u])}

            non_sep_sp = are_species_sep(components)
            algo = 'GN'


            #If not all the species are separated try kerningan_lin
            if non_sep_sp >= 1:

                for _ in range(5): #5 random starts

                    #FIXME we could try to get a deterministic result with the KL algo:
                    #from tests it seems specifying a random seed to this function does not guarantee
                    #determinism
                    bisect = nx.algorithms.community.kernighan_lin_bisection(graph, partition=None,
                                                                             max_iter=10,
                                                                             weight='weight')

                    cut_edges = {(u, v) for u in bisect[0]\
                                 for v in bisect[1].intersection(graph.adj[u])}

                    non_sep = are_species_sep(bisect)

                    #if better separation keep KL
                    if non_sep < non_sep_sp:
                        cut_edges_best = cut_edges
                        components = bisect
                        non_sep_sp = non_sep
                        algo = 'KL'

    return components, cut_edges_best, algo, (scores_cut, scores_uncut)


def worker_cut_graph(family, fam, res, spectral=False, k=0, verbose=False):
    """
    Worker for parallel graph cutting. Collapses tandem duplicates, detects the two communities in
    the graph and store results and statistics about the cuts in the `res` dictionary.

    Args:
        family (networkx.Graph): Orthology graph of the gene family
        fam (str): Unique id of the gene family
        res (dict): Dictionary storing the results, shared between processes
        spectral (bool, optional): Use spectral clustering instead of default Girvan-Newman
                                   (faster)
        k (int, optional): unique id for the cut graph
        verbose (bool, optional): print progress

    Returns:
        bool: True if no Exception was raised.
    """

    try:
        #some nodes are place-holder for a loss of duplicates, called 'None_nb_species'
        #where `nb` is the index of the family.
        n_species = len({species_name(gene) for gene in family\
                         if gene != 'None_'+species_name(gene)})
        n_genes = len([gene for gene in family if gene != 'None_'+species_name(gene)])


        #Filter multigenic families
        if n_genes < n_species*3:

            #collapse tandem duplicates
            collapse_nodes(family)

            #cut graphs
            orthogroups, cuts, algo, scores = min_cut(family, spectral)

            #ignore graphs with two genes
            if algo == 'Too few genes':
                return True

            orthogroups = list(orthogroups)

            #store orthogroups
            f_ortho = []
            for sublist in orthogroups:

                new = []

                #filter loss and extend the groups by adding collapsed node
                for i in sublist:

                    if i != 'None_'+species_name(i):

                        if 'contraction' in family.node[i]:

                            new.extend([gene for gene in family.node[i]['contraction']])

                        new.append(i)

                f_ortho.append(new)

            res[fam] = (f_ortho, algo, len(cuts), scores)

        else:
            res[fam] = (None, 'Filtered_Multigenic', 'NAN', 'None')

        if verbose:
            sys.stderr.write(f"Cut graph {k} \n")
            sys.stderr.flush()

        return True

    except Exception:

        traceback.print_exc()
        raise

def print_out_stats(stats_dict, wgd=''):

    """
    Prints to stdout some statistics about community detection in graphs.

    Args:
        stats_dict (dict): a dict counting the number of graphs that were either multigenic and
                           discarded or processed and cut with different algo
        wgd (str, optional): the wgd for which the synteny graphs were computed
    """

    if stats_dict:

        multi = stats_dict.get("Filtered_Multigenic", 0)
        gn = stats_dict.get("GN", 0)
        kl = stats_dict.get("KL", 0)
        cl = stats_dict.get("No", 0)
        tot = gn + kl + cl

        if tot:

            print('\n')
            print("-----------------------COMMUNITY DETECTION IN GRAPHS------------------------")
            print(" Whole-genome duplication: {}".format(wgd))
            print("\n")
            print(" {} total orthology graphs effectively built".format(tot+multi))
            print(" (Families can fail to produce graphs if gene members cannot be threaded "
                  "\n in the synteny analysis)")
            print("\n")
            print(" {} processed graphs out of {} ({} discarded multigenic families)"\
                  .format(tot, tot+multi, multi))

            clp = round((cl/tot)*100, 2)
            gnp = round((gn/tot)*100, 2)
            klp = round((kl/tot)*100, 2)

            print(" {} ({} %) graphs were two separated cliques".format(cl, clp))

            print(" {} ({} %) graphs cut with the Girvan-Newman algo.".format(gn, gnp))
            print(" {} ({} %) graphs cut with the Kerningan-Lin algo.".format(kl, klp))
            print("----------------------------------------------------------------------------")
            print("\n")

        else:
            print('\n')
            print("-----------------------COMMUNITY DETECTION IN GRAPHS------------------------")
            print(" Whole-genome duplication: {}".format(wgd))
            print("\n")
            print(" {} processed graphs out of {} ({} discarded multigenic families)"\
                  .format(tot, tot+multi, multi))
            print("----------------------------------------------------------------------------")
            print("\n")


if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required

    PARSER.add_argument('-i', '--input', type=str, help='Single file with all orthologies\
                        prediction, sorted and gzipped.', required=True)

    #Optional

    PARSER.add_argument('-o', '--output', type=str, help='Output file', required=False,
                        default="out")

    PARSER.add_argument('-w', '--weight', type=str, help='Should weights be used in the orthology\
                        graphs?', required=False, default='n')

    PARSER.add_argument('-n', '--ncores', type=int, help='Number of threads for the parallel\
                        cutting of graphs.', required=False, default=1)

    PARSER.add_argument('-s', '--summary', type=str, help='File to write a summary on cuts',
                        required=False, default='Summary')

    PARSER.add_argument('-ignSg', '--ignoreSingleGeneC', type=str,
                        help='Whether to ignore single gene communities',
                        required=False, default='y')

    PARSER.add_argument('-wgd', '--wgd_tag', type=str,
                        help='Informations on the WGD to print it out along with statistics',
                        required=False, default='')

    PARSER.add_argument('--spectral', action='store_true', help="Force community detection"
                        "with spectral clustering, for computational efficiency")

    PARSER.add_argument('--verbose', action='store_true')

    ARGS = vars(PARSER.parse_args())

    for ARG in ['ignoreSingleGeneC', 'weight']:
        assert ARGS[ARG] in ['y', 'n'], '{} should be y or n'.format(ARG)

        if ARGS[ARG] == 'y':
            ARGS[ARG] = True
        else:
            ARGS[ARG] = False

    #Lazy load each graph and find the 2 post-duplication communities in each family.
    sys.stderr.write('Cutting the orthology graphs...\n')

    #dictionary shared between processes
    MANAGER = multiprocessing.Manager()
    RES = MANAGER.dict()

    POOL = multiprocessing.Pool(ARGS["ncores"])
    ASYNC_RESULTS = []
    with gzip.open(ARGS["input"], 'rt') as infile:
        k = 0
        for FAM, FAM_ID in lazy_load_pairwise_file(infile, use_weights=ARGS['weight']):
            JOB = POOL.apply_async(worker_cut_graph, args=(FAM, FAM_ID, RES, ARGS["spectral"], k,
                                                           ARGS["verbose"]))
            ASYNC_RESULTS += [JOB]
            k += 1

    POOL.close()
    POOL.join()

    for RESULT in ASYNC_RESULTS:
        if not RESULT.get():
            sys.stderr.write("An error occured in a child process\n")
            sys.exit(1)

    #Write results to file
    # ALL_SCORES = []
    STATS = {}
    sys.stderr.write('Writting graphs orthogroups...\n')

    with open(ARGS['output'], 'w') as OUTFILE, open(ARGS['summary'], 'w') as SUMMARY:

        SUMMARY.write('Family\tCommunity detection algo.\tNumber of removed edges\n')

        for FAM in RES.keys():

            ORTHOGROUPS, ALGO, CUTS, SCORES = RES[FAM]

            if ORTHOGROUPS or ALGO == 'Filtered_Multigenic':

                SUMMARY.write(FAM+'\t'+ALGO+'\t'+str(CUTS)+'\n')
                STATS[ALGO] = STATS.get(ALGO, 0) + 1

            #write first community (largest)
            if ORTHOGROUPS:
                ORTHOGROUPS.sort(key=len, reverse=True)
                if (len(ORTHOGROUPS[0]) > 1 and ARGS["ignoreSingleGeneC"]) or\
                   (len(ORTHOGROUPS[0]) >= 1 and not ARGS["ignoreSingleGeneC"]):

                    ORTHO_A = ('\t').join([FAM+'a'] + list(ORTHOGROUPS[0])) #arbitrary a
                    OUTFILE.write(ORTHO_A+'\n')
                    # ALL_SCORES.append(SCORES)

                #write second community (if more than two in total, write the second largest)
                if len(ORTHOGROUPS) > 1 and len(ORTHOGROUPS[1]) >= 2 and ARGS["ignoreSingleGeneC"]:
                    ORTHO_B = ('\t').join([FAM+'b'] + list(ORTHOGROUPS[1])) #arbitrary b
                    OUTFILE.write(ORTHO_B+'\n')

                elif len(ORTHOGROUPS) > 1 and len(ORTHOGROUPS[1]) >= 1\
                                          and not ARGS["ignoreSingleGeneC"]:
                    ORTHO_B = ('\t').join([FAM+'b'] + list(ORTHOGROUPS[1])) #arbitrary b
                    OUTFILE.write(ORTHO_B+'\n')

    if ARGS["wgd_tag"]:
        WGD, OUTGR = ARGS["wgd_tag"].split(',')
        WGD_TAG = WGD + ' (outgroup '+OUTGR+')'
        print_out_stats(STATS, wgd=WGD_TAG)
