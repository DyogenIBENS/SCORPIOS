import treeCl
import argparse
import os

from collections import defaultdict
import csv
import functools
import itertools
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


def load_convert_distmat(input_file, output_file):
    # name_to_id associates a name with an integer 0, 1, ...
    name_to_id = defaultdict(functools.partial(next, itertools.count()))

    with open(input_file) as infile:
        reader = csv.reader(infile)

        #store row and column names in dict
        for name_a, name_b, dist in reader:
            idx_a = name_to_id[name_a]
            idx_b = name_to_id[name_b]

        # make the distance matrix
        n_elem = len(name_to_id)
        dists = np.zeros((n_elem, n_elem))

        # go back to the start of the file and load the data
        infile.seek(0)
        for name_a, name_b, dist in reader:
            idx_a = name_to_id[name_a]
            idx_b = name_to_id[name_b]

            if dist == np.nan or dist == "nan":
                dist = 0.0

            dists[(idx_a, idx_b)] = dist
            dists[(idx_b, idx_a)] = dist


    id_to_name = dict((identifier, name) for name, identifier in name_to_id.iteritems())

    with open(output_file, 'w') as fw:
        fw.write(','.join([id_to_name[k] for k in range(len(id_to_name))])+'\n')
        for i in range(len(id_to_name)):
            name = id_to_name[i]

            d = ','.join([str(dists[(i, k)]) for k in range(len(id_to_name))])
            fw.write(name+','+d+'\n')





if __name__ == '__main__':
    

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-d', '--distmat', required=True)

    PARSER.add_argument('-od', '--outdistmat', required=True)

    PARSER.add_argument('-om', '--outmds', required=True)

    PARSER.add_argument('-o', '--output', required=True)

    PARSER.add_argument('-k', '--k', required=False, default=3, type=int)

    ARGS = vars(PARSER.parse_args())

    OUT, EXT = os.path.splitext(ARGS["distmat"])
    OUTDISTMAT = OUT+"_converted"+EXT

    load_convert_distmat(ARGS["distmat"], OUTDISTMAT)


    OUTFIGDM = ARGS["outdistmat"]
    OUTFIGMDS = ARGS["outmds"]
    OUTCLUST = ARGS["output"]

    # for i in range(2, 5):
    dm = treeCl.DistanceMatrix.from_csv(OUTDISTMAT)
    # print(dm)
    spclust = treeCl.Spectral(dm)
    # print(spclust)
    partitions = spclust.cluster(ARGS["k"])
    treeCl.plotter.heatmap(dm, partition=partitions)
    plt.savefig(OUTFIGDM, dpi=100)
    # plt.show()
    coord = spclust.spectral_embedding(3)
    # print(treeCl.Evaluation(dm).silhouette(partitions))
    treeCl.plotter.plot_embedding(coord, partition=partitions)
    plt.savefig(OUTFIGMDS, dpi=100)
    # plt.show()

    plt.close("all")
    loci = dm.get_names()
    clusters = partitions.partition_vector
    with open(OUTCLUST, 'w') as fw:
        for tree, clust in zip(loci, clusters):
            fw.write(tree+'\t'+str(clust)+'\n')

    #grep nan full_matrix_RF_4000sample.csv | grep -Eo ENSELUG[0-9]+ | sort | uniq -c | awk -v limit=10 '$1 > limit{print $2}' > problematic_treees
