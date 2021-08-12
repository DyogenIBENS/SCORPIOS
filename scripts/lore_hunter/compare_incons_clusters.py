import sys
import argparse

from .medoid_topologies import load_clust

if __name__ == '__main__':
        # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-c', '--clusters', required=True)

    PARSER.add_argument('-o', '--output', required=True)

    PARSER.add_argument('-i', '--incons', required=True)

    ARGS = vars(PARSER.parse_args())

    CLUSTERS = load_clust(ARGS["clusters"])
    INCONS = load_clust(ARGS["incons"], ["Inconsistent"])

    ALL_CLUSTERS = {tree for i in CLUSTERS for tree in CLUSTERS[i]}
    TOT = len(set(INCONS["Inconsistent"]).intersection(ALL_CLUSTERS))

    with open(ARGS["output"], 'w') as OUT:
        OUT.write(f"{TOT} inconsistent trees sampled for clustering:\n")
        sys.stdout.write(f"{TOT} inconsistent trees sampled for clustering:\n")
        for i in CLUSTERS:

            incons = len([tree for tree in CLUSTERS[i] if tree in INCONS["Inconsistent"]])
            lg = len(CLUSTERS[i])
            OUT.write(f"{incons} inconsistent trees in cluster {i} (cluster of size={lg} trees)\n")
            sys.stdout.write(f"{incons} inconsistent trees in cluster {i} (cluster of size={lg} trees)\n")

