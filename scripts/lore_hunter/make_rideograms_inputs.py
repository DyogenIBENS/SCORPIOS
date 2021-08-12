"""

"""
import sys

import roman
import argparse

from scripts.synteny.mygenome import Genome

def make_karyo(genesfile, output, fomt='bed'):
    """
    """
    genome = Genome(genesfile, fomt)
    dgenes = genome.genes_list
    karyo = []
    for chrom in dgenes:

        if isinstance(chrom, str) and "HiC" in chrom:
            new_chrom = chrom.replace("HiC_scaffold_", "")

        elif isinstance(chrom, str) and "chr" in chrom:
            new_chrom = chrom.replace("chr", "")

        elif isinstance(chrom, str) and "group" in chrom:
            new_chrom = chrom.replace("group", "")

        elif isinstance(chrom, str) and "LG" in chrom:
            new_chrom = chrom.replace("group", "")

        else:
            new_chrom = chrom

        try:
            new_chrom = roman.fromRoman(new_chrom)
        except:
            pass

        try:
            karyo.append((int(new_chrom), len(dgenes[chrom]), chrom))
        except:
            sys.stderr.write(f"Warning: chromosome {chrom} could not be converted to an integer, it won't be drawn.")
            continue

    karyo = sorted(karyo)

    with open(output, 'w') as outfile:
        outfile.write("Chr\tStart\tEnd\n")
        
        for ch in karyo:
            chrom, lg, str_chrom = ch
            outfile.write(f"{chrom}\t0\t{lg+1}\n")
        
    return genome, karyo


def load_features(genome, features_file, to_load=None):
    feat = {}
    sp_genes = {g.names[0] for g in genome}
    with open(features_file, 'r') as infile:

        for line in infile:

            _, descendants, classif = line.strip().split('\t')


            if to_load is None or classif in to_load:

                try:
                    classif = int(classif)
                except:
                    if classif == "Inconsistent":
                        classif = 1
                    else:
                        pass

                genes = set(descendants.split()).intersection(sp_genes)


                if genes:
                    genes = list(genes)
                    for g in genes:
                        feat[g] = classif
    return feat

def features_to_ide(genome, features_file, karyo, output):

    feat = load_features(genome, features_file)

    dgenes = genome.genes_list

    with open(output, 'w') as out:
        out.write("Chr\tStart\tEnd\tValue\n")

        for k in karyo:
            chrom_int, lg, chrom = k
            prev = False

            for i, gene in enumerate(dgenes[chrom]):
                name = gene.names[0]
                to_write = False

                if name in feat:
                    # if feat[name] == "Inconsistent":
                    classif = feat[name]
                    to_write = True

                if to_write and not prev:
                    start = i
                    stop = i + 1
                    prev = feat[name]

                elif to_write and prev == classif:
                    stop = i + 1

                elif (not to_write and prev) or (to_write and prev != classif):
                    out.write(f"{chrom_int}\t{start}\t{stop}\t{prev}\n")
                    prev = False
                    if name in feat:
                        prev = classif
                        start = i
                        stop = i + 1

            if prev:
                out.write(f"{chrom_int}\t{start}\t{stop}\t{prev}\n")
                prev = False

if __name__ == '__main__':

    # Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('-i', '--input_incons', help='', required=True)

    PARSER.add_argument('-g', '--genesfile', help='', required=True)

    PARSER.add_argument('-f', '--format', help='', required=False, default="bed")

    PARSER.add_argument('-k', '--outfile_karyo', help='Output file rideogram karyotype.', required=False, default="out")

    PARSER.add_argument('-o', '--outfile_features', help='Output file rideogram overlaid data.', required=False, default="out")

    ARGS = vars(PARSER.parse_args())

    GENOME, KARYO = make_karyo(ARGS["genesfile"], ARGS["outfile_karyo"], fomt=ARGS["format"])

    features_to_ide(GENOME, ARGS["input_incons"], KARYO, ARGS["outfile_features"])