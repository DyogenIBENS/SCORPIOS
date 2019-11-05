#!/usr/bin/env python

"""
    Script to find all orthology relationships between a group of WGD duplicated species and a non
    duplicated outgroup. These ortholog groups define gene families.

    Example:
        $ python -m scripts.synteny.duplicated_families -t forest_v89.nhx -n Lepisosteus.oculatus
                                                        -d Clupeocephala -s Species_tree_v89.nwk
                                                        -g genes89/genesST.%s.list.bed [-o out]
                                                        [-ow anc1,anc2] [-u ufile]
"""

import sys
import argparse
import itertools

from ete3 import Tree

from scripts.trees import speciestree as spt, utilities as ut, orthologs as org
from . import mygenome
from . import utilities as syu


def tag_duplicated_species(leaves, duplicated):

    """
    Adds a tag to genes of duplicated species in an ete3 Tree instance, in-place.

    Args:
        leaves (list of ete3 TreeNode): Leaves of the tree.
        duplicated (list of str): List of the names of all duplicated species
    """

    for leaf in leaves:

        if leaf.S in duplicated:
            leaf.add_features(duplicated='Y')

        else:
            leaf.add_features(duplicated='N')


def get_genes_positions(genes, species, dict_genes):

    """
    Gets genomic position of given `genes` of a species.

    Args:
        sp (str): Input species name.
        genes (list of str): List of the genes to search.
        dict_genes (dict of str:GeneSpeciesPosition tuples): Genes location.

    Returns:
        list: Genes and their position as a list of GeneSpeciesPosition tuples
    """

    ortho_genes = []
    for gene in genes:
        if gene in dict_genes[species]:
            ortho_genes.append(dict_genes[species][gene])

        else:
            sys.stderr.write("Warning: {} is not in the genes coordinate file \n".format(gene))
            ortho_genes.append(syu.GeneSpeciesPosition(gene, '', ''))

    return ortho_genes


def orthologies_with_outgroup(forest, duplicated_sp, outgroup, dict_genes, out):

    """
    Browses a gene tree forest and searches for orthologs with the outgroup.
    Writes genes without phylogenetic orthologs to a file.

    Args:
        forest (str): Name of the gene trees forest file
        duplicated_sp (list of str): List of all duplicated species for the considered WGD
        outgroup (str): Non-duplicated outgroup
        dict_genes (dict of GeneSpeciesPosition tuples): All gene positions for each species
        out (str): Output file to write genes without phylogenetic orthologs

    Returns:
        dict: Orthologs of outgroup genes in each duplicated species

    """

    ortho = {e: {} for e in duplicated_sp}

    with open(out, 'w') as outfile, open(forest, 'r') as infile:

        sys.stderr.write("Browsing gene trees for orthologies with the outgroup...\n")

        for tree in ut.read_multiple_objects(infile):

            #load tree
            tree = Tree(tree.strip(), format=1)
            node2leaves = tree.get_cached_content()
            leaves = [i for i in tree.get_leaves()]

            #add a tag to genes of duplicated species
            tag_duplicated_species(leaves, duplicated_sp)

            #find all clades with only genes of duplicated species
            subtrees = tree.get_monophyletic(values=["Y"], target_attr="duplicated")

            #find all outgroup genes
            outgroup_genes = [i for i in leaves if i.S == outgroup]

            #search for an ortholog gene in the outgroup for all clades of teleost genes
            for subtree in subtrees:

                seen = {}
                subtree_leaves = subtree.get_leaves()
                found = False

                #browse all outgroup genes
                for j in outgroup_genes:

                    #find the node that splits the outgroup gene and duplicated species genes
                    lca = tree.get_common_ancestor(subtree, j)
                    topo_distance = len(node2leaves[lca])

                    # if it is a speciation or dubious duplication node --> speciation
                    if org.is_speciation(lca):
                        branch_distance = tree.get_distance(subtree, j)
                        if subtree not in seen:
                            seen[subtree] = []
                        seen[subtree].append((topo_distance, branch_distance, j))
                        found = True
                        #break

                # if no 'true' ortholog
                # check if all descendants include only outgroup + duplicated species
                if not found:
                    for j in outgroup_genes:
                        lca = tree.get_common_ancestor(subtree, j)

                        for gene in lca.get_leaves():
                            if gene.duplicated != "Y" and gene.S != outgroup:
                                break

                        #if no break, it means all descendants are outgroup or dup.
                        else:
                            topo_distance = len(node2leaves[lca])
                            branch_distance = tree.get_distance(subtree, j)
                            seen[subtree] = seen.get(subtree, [])
                            seen[subtree].append((topo_distance, branch_distance, j))


                # if an ortholog was found, add it to the orthology dict
                if seen:
                    seen[subtree].sort(key=lambda x: (x[0], x[1]))
                    outgroup_gene = seen[subtree][0]

                    outgroup_gene = outgroup_gene[2].name
                    for species in duplicated_sp:
                        genes = [i.name for i in subtree_leaves if i.S == species]
                        genes = get_genes_positions(genes, species, dict_genes)

                        ortho[species][outgroup_gene] = ortho[species].get(outgroup_gene, [])
                        ortho[species][outgroup_gene] += genes


                # if no ortholog found
                # write genes without ortholog along with all outgroup genes in tree
                # (potential candidate for orthology)
                elif any(i.name in dict_genes[outgroup] for i in outgroup_genes):

                    #genes without orthologs
                    missed_genes = []
                    for species in duplicated_sp:
                        genes = [i.name for i in subtree_leaves if i.S == species]
                        genes = get_genes_positions(genes, species, dict_genes)
                        missed_genes += [g.name+'_'+species.replace(' ', '.')+\
                                         '|'+str(g.chromosome)+\
                                         '|'+str(g.index) for g in genes]

                    if missed_genes:
                        outfile.write(' '.join(missed_genes)+'\t')

                        #candidate orthologs in the outgroup
                        outgr_genes = [i.name for i in outgroup_genes]

                        in_paralogs = []
                        for pair in itertools.combinations(outgr_genes, 2):
                            if tree.get_distance(pair[0], pair[1], topology_only=True) == 1:
                                in_paralogs.append(pair[0]+'|'+pair[1])

                        outgr_write = []
                        genome = dict_genes[outgroup]
                        for gene in outgr_genes:
                            if gene in genome:

                                lca = tree.get_common_ancestor(subtree, gene)
                                branch_distance = tree.get_distance(subtree, gene)
                                topo_distance = len(node2leaves[lca])

                                outgr_write.append(str(gene)+'|'+str(genome[gene].chromosome)+'|'+\
                                                   str(genome[gene].index)+'|'+str(topo_distance)+\
                                                   '|'+str(branch_distance))

                        outfile.write(' '.join(outgr_write)+'\t'+' '.join(in_paralogs)+'\n')

    sys.stderr.write("Phylogenetic orthologies with the outgroup OK\n")

    return ortho


def write_orthologs(orthos, dicgenomes, dict_genes, outgroup, duplicated_sp, out, min_length=20):

    """
    Writes to a file gene orthologies between the non-duplicated species and all duplicated
    species (orthologytable), with all gene names and gene positions. All these gene families are
    ordered along the outgroup genome in the output.

    Args:
        orthos (dict of str:str:GeneSpeciesPosition tuples): orthologs of outgroup genes in each
                                                             duplicated species.
        dicgenomes (dict of str:mygenome.Genome): Genomes.
        dict_genes (dict of str:GeneSpeciesPosition tuples): Genes location.
        outgroup (str): non-duplicated outgroup
        duplicated_sp (list of str): List of duplicated species to include in the results.
        out (str): Output file name for genes without orthologs
        min_length (int, optional): Minimum length for a chromosome in the outgroup, gene families
                                    mapping to smaller chromosomes won't be included.
    """

    sys.stderr.write("Writing phylogenetic orthologies with the outgroup...\n")

    stats = {}

    with open(out, 'w') as outfile:

        outfile.write('\t'.join([outgroup, "", "", ""] + duplicated_sp)+'\n')

        genome = dicgenomes[outgroup]

        # Chromosomes of the non-duplicated species
        for chrom in genome.chr_list[mygenome.ContigType.Chromosome] +\
                     genome.chr_list[mygenome.ContigType.Scaffold]:

            # Do not write ortholog for very small scaffolds
            if len(genome.genes_list[chrom]) < min_length:
                continue

            #Genes on the chromosome of the non-duplicated species
            for gene in genome.genes_list[chrom]:

                name = gene.names[0]

                (_, _, i) = dict_genes[outgroup][name]

                all_orthologs = []
                found_one = False

                for spd in duplicated_sp:

                    stats[spd] = stats.get(spd, 0)

                    orthologs = [""]

                    # We need orthologs in duplicated species
                    if name in orthos[spd]:
                        orthologs = orthos[spd][name]
                        found_one = True
                        stats[spd] += len(orthologs)
                        orthologs = [('/'.join(["{}|{}-{}".format(getattr(o, "name"),\
                                      getattr(o, "chromosome"), getattr(o, "index"))\
                                      for o in orthologs]))]
                    all_orthologs += orthologs

                #write if at least one ortholog in a duplicated species
                if found_one:
                    stats["families"] = stats.get("families", 0)
                    stats["families"] += 1
                    all_orthologs = ["\t".join(all_orthologs)]
                    outfile.write('\t'.join([str(chrom), str(i), name, ""] + all_orthologs)+'\n')

    return stats


def print_out_stats(stats_dict, wgd=''):

    """
    Prints to stdout some statistics on the families in the phylogenetic Orthology Table.

    Args:
        stats_dict (dict): a dict counting number of families and genes in the families
        wgd (str, optional): the wgd for which the Orhtology Table was built
    """

    empty_table = True

    if stats_dict:

        fam = stats_dict.get('families', 0)

        if fam:

            empty_table = False

            print('\n')
            print("------------------------PHYLOGENETIC ORTHOLOGY TABLE------------------------")
            print(" Whole-genome duplication: {}".format(wgd))
            print("\n")
            print(" {} total families in the phylogenetic orthology table".format(fam))
            print("\n")
            for spec in stats_dict:
                if spec != 'families':
                    genes = stats_dict[spec]
                    print(" {} {} genes in the table".format(genes, spec))
            print("----------------------------------------------------------------------------")
            print("\n")

    if empty_table:
        print('\n')
        print("------------------------PHYLOGENETIC ORTHOLOGY TABLE------------------------")
        print(" Whole-genome duplication: {}".format(wgd))
        print("\n")
        print(" 0 total families in the phylogenetic orthology table")
        print("----------------------------------------------------------------------------")
        print("\n")


if __name__ == '__main__':


    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-t', '--treesFile', help='Forest of trees in .nhx, with species,\
                         duplication/speciation nodes.',\
                         required=True)

    PARSER.add_argument('-n', '--outgroup', help='Non-duplicated outgroup.', required=True)

    PARSER.add_argument('-d', '--dupSp', help='Name of the ancestor of all duplicated species.',\
                         required=True)

    PARSER.add_argument('-s', '--speciesTree', help='Species tree (newick), with ancestor names.',\
                         required=True)

    PARSER.add_argument('-g', '--genesFile', help='Genomic coordinates of genes (.bed) or (.bz2)',\
                         required=True)

    PARSER.add_argument('-o', '--out', help='Output file', required=False, default="out")

    PARSER.add_argument('-ow', '--dupSp2', help='ancestor(s) with other WGDs', required=False)

    PARSER.add_argument('-u', '--uncertainFamilies', help='File to write genes without\
                         phylogenetic orthologs.', required=False, default="Uncertain_genes")

    ARGS = vars(PARSER.parse_args())

    #Study species
    OTHER_WGD = ""
    if 'dupSp2' in ARGS and ARGS["dupSp2"]:
        OTHER_WGD = ARGS["dupSp2"]
    DUP_SPECIES = spt.get_species(ARGS["speciesTree"], ARGS["dupSp"], OTHER_WGD)
    OUTGROUP = ARGS["outgroup"]

    #Get genes positions
    DICTGENOME, DICTGENES = {}, {}
    SPECIES = DUP_SPECIES + [OUTGROUP]
    for SP in SPECIES:
        GENOME = mygenome.Genome(ARGS["genesFile"] % SP)
        if SP == OUTGROUP:
            DICTGENOME[SP] = GENOME
        DICTGENES[SP] = {}
        for (NAME, (CHROM, INDEX)) in GENOME.dict_genes.items():
            DICTGENES[SP][NAME] = syu.GeneSpeciesPosition(NAME, CHROM, INDEX)


    #Extract orthologies from gene trees
    ORTHOS = orthologies_with_outgroup(ARGS["treesFile"], DUP_SPECIES, OUTGROUP,\
                                       DICTGENES, ARGS["uncertainFamilies"])

    #Write orthology file (orthology mapping of duplicated species on the outgroup)
    STATS = write_orthologs(ORTHOS, DICTGENOME, DICTGENES, OUTGROUP, DUP_SPECIES, ARGS["out"])

    print_out_stats(STATS, wgd=ARGS["dupSp"]+' (Outgroup '+ARGS["outgroup"]+')')
