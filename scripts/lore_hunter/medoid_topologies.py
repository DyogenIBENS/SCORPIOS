import argparse
import os

from ete3 import Tree
import pandas as pd

def load_matrix(input_file):
    matrix = pd.read_csv(input_file, index_col=0)
    return matrix

def get_n_medoids(matrix, clusters, n=5):
    res = {}

    for clust in clusters:
        tmp_mat = matrix.loc[matrix.index.intersection(clusters[clust])]
        tmp_mat["sum"] = tmp_mat.sum(axis=1)
        medoids = tmp_mat.nsmallest(n, ["sum"], keep='first').index.tolist()
        res[clust] = medoids
    return res

def load_clust(input_file, load_only=None):
    with open(input_file, 'r') as infile:
        clust = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in infile}

    res = {}
    for tree in clust:
        clust_id = clust[tree]
        if load_only is None or clust_id in load_only:
            res[clust_id] = res.get(clust_id, [])
            res[clust_id].append(tree)
    return res

def write_medoids(dict_medoids, trees_path, outfolder, load_only=None):
    os.makedirs(f"{outfolder}", exist_ok=True)
    for clust in dict_medoids:
        for i, medoid in enumerate(dict_medoids[clust]):
            tree = Tree(trees_path.format(medoid))
            if not load_only: #dirty hack, to improve
                tree.write(outfile=f"{outfolder}/cluster_{clust}_medoid_{i+1}.nhx", format=1, features=["S", "D"])
            else:
                tree.write(outfile=f"{outfolder}/medoid_{i+1}.nhx", format=1, features=["S", "D"])


if __name__ == '__main__':
    

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-d', '--distmat', required=True)

    PARSER.add_argument('-c', '--clusters', required=True)

    PARSER.add_argument('-o', '--outdir', required=True)

    PARSER.add_argument('-t', '--trees_path', required=True)

    PARSER.add_argument('-n', '--n', required=False, type=int, default=5)

    PARSER.add_argument('-r', '--restrict', required=False, nargs="*", default=None)

    ARGS = vars(PARSER.parse_args())

    MAT = load_matrix(ARGS["distmat"])
    CLUSTERS = load_clust(ARGS["clusters"], ARGS["restrict"])
    MED = get_n_medoids(MAT, CLUSTERS, ARGS["n"])
    write_medoids(MED, ARGS["trees_path"], ARGS["outdir"], ARGS["restrict"])
