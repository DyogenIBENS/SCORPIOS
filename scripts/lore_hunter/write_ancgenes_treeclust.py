"""

"""

import os.path
import argparse
from ete3 import Tree


def write_ancgenes(clustered_genes, treedir, out_ancgenes, clusters_to_load = None):

    """

    """

    k = 0

    with open(out_ancgenes, 'w') as outfile:

        for gene in clustered_genes:

            cluster = clustered_genes[gene]

            if clusters_to_load is not None and cluster not in clusters_to_load: #!= "Inconsistent"
                continue

            treefile = treedir +  '/' + gene + '.nhx'

            if not os.path.exists(treefile):
                treefile = treedir +  '/' + gene + '.nh'

            if not os.path.exists(treefile):

                treefile = treedir +  '/C_' + gene + '.nh'

            if not os.path.exists(treefile):
                treefile = treedir + "/" + gene + "_final.nhx"

            assert os.path.exists(treefile), f"The file {treefile} does not exist"

            tree = Tree(treefile)

            leaves = {'_'.join(i.name.split('_')[:-1]) for i in tree.get_leaves()}

            if leaves == {''}:
                leaves = {i.name for i in tree.get_leaves()}
            
            descendants = sorted(list(leaves))

            anc = 'Name_'+str(k)

            if clusters_to_load is not None:
                cluster = str(clusters_to_load.index(cluster))
                    
            outfile.write(anc+'\t'+ ' '.join(descendants)+'\t'+cluster+'\n')

            k += 1


def load_gene_list(input_summary, input_acc=None):

    """

    """

    with open(input_summary, 'r') as infile:
        genes = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in infile}

    if input_acc is not None:

        with open(input_acc, 'r') as infile:
            acc = {line.strip().split('\t')[0] for line in infile}

        for g in genes:
            if genes[g] == "Inconsistent" and g in acc:
                genes[g] = "Consistent"

    return genes


def write_summary(summary_dict, output_file):
    """
    """
    with open(output_file, 'w') as out:
        for key in summary_dict:
            out.write('\t'.join([key, summary_dict[key]])+'\n')


if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--treesdir', help='', required=True)

    PARSER.add_argument('-c', '--clusters', help='', required=True)

    PARSER.add_argument('-a', '--accepted', help='SCORPiOs accepted corrections',
                        required=False, default=None)

    PARSER.add_argument('-o', '--outfile', help='Output file', required=False, default="out")

    PARSER.add_argument('--summary_only', help='Only write outgroup gene name + tree consistency.',
                        action="store_true")

    PARSER.add_argument('-r', '--restrict_to', required=False, default=None, nargs='*')

    ARGS = vars(PARSER.parse_args())

    CLUSTERS = load_gene_list(ARGS["clusters"], ARGS["accepted"])

    if ARGS["summary_only"]:
        write_summary(CLUSTERS, ARGS["outfile"])

    else:
        write_ancgenes(CLUSTERS, ARGS["treesdir"], ARGS["outfile"], ARGS["restrict_to"])
