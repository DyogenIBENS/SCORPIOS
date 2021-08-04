"""

"""

import os.path
import argparse
from ete3 import Tree


def write_ancgenes(clustered_genes, treedir, out_ancgenes):

    """

    """

    k = 0

    with open(out_ancgenes, 'w') as outfile:

        for gene in clustered_genes:

            cluster = clustered_genes[gene]

            if cluster != "Inconsistent":
                continue

            treefile = treedir +  '/' + gene + '.nhx'

            if not os.path.exists(treefile):

                treefile = treedir +  '/C_' + gene + '.nh'

            assert os.path.exists(treefile), f"The file {treefile} does not exist"

            tree = Tree(treefile)

            leaves = {'_'.join(i.name.split('_')[:-1]) for i in tree.get_leaves()}
            
            descendants = sorted(list(leaves))


            anc = 'Name_'+str(k)
                    
            outfile.write(anc+'\t'+ ' '.join(descendants)+'\t'+cluster+'\n')

            k += 1


def load_gene_list(input_summary, input_acc):

    """

    """

    with open(input_summary, 'r') as infile:
        genes = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in infile}

    with open(input_acc, 'r') as infile:
        acc = {line.strip().split('\t')[0] for line in infile}

    for g in genes:
        if genes[g] == "Inconsistent" and g in acc:
            genes[g] = "Consistent"

    return genes


if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--treesdir', help='', required=True)

    PARSER.add_argument('-c', '--clusters', help='', required=True)

    PARSER.add_argument('-a', '--accepted', help='SCORPiOs accepted corrections',
                        required=True)


    PARSER.add_argument('-o', '--outfile', help='Output file', required=False, default="out")


    ARGS = vars(PARSER.parse_args())

    CLUSTERS = load_gene_list(ARGS["clusters"], ARGS["accepted"])

    write_ancgenes(CLUSTERS, ARGS["treesdir"], ARGS["outfile"])
