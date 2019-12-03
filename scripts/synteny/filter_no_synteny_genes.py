"""
    This script identifies genes in the orthology table that never, in any of their sliding
    windows, have genes on the same chromosome in the orthology table.
    A new orthology table is written as output, where genomic posistion of these genes is omitted,
    which forces SCORPiOs other scripts to not use them in the synteny analysis.

    Example:
        $ python -m scripts.synteny.filter_no_synteny_genes -i OrthoTable.txt -chr Chr_outgr_file
                                                            [-o out] [-w 15]
"""


import argparse
import collections

import numpy as np

from . import syntenycompare as synt
from . import utilities as ut

if __name__ == '__main__':

    ## Arguments
    PARSER = argparse.ArgumentParser(description=__doc__,\
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Required

    PARSER.add_argument('-i', '--input', type=str, help='Orthology Table', required=True)

    PARSER.add_argument('-chr', '--chr_outgr', type=str, help='Outgroup chromosome to consider.',
                        required=True)

    # Optional

    PARSER.add_argument('-o', '--out', type=str, help='Result file', required=False, default="out")

    PARSER.add_argument('-w', '--windowSize', type=int, help='Size of the sliding window',
                        required=False, default=15)

    ARGS = vars(PARSER.parse_args())

    ORTHOTABLE = ARGS["input"]

    with open(ORTHOTABLE, 'r') as infile:
        FIRST_LINE = infile.readline()
        SPECIES = [i for i in FIRST_LINE.strip().split('\t')[1:] if i]
        OUTGR = FIRST_LINE.strip().split('\t')[0]

    CHROMOSOMES = []

    with open(ARGS["chr_outgr"], 'r') as infile:
        CHROMOSOMES += [line.strip() for line in infile]


    DGENES = {}
    ALL_ORTHOS = collections.defaultdict(dict)

    for chrom in CHROMOSOMES:
        for spec in SPECIES:

            if spec not in DGENES:
                DGENES[spec] = {"synteny":set(), "nosynteny":set()}

            table_entries_sp = ut.complete_load_orthotable(ORTHOTABLE, chrom, spec)

            ALL_ORTHOS[chrom][spec] = ut.complete_load_orthotable(ORTHOTABLE, chrom, spec,
                                                                  load_no_position_genes=True)

            start = 0
            stop = len(table_entries_sp)

            #check that regions is at least as long as windowSize
            if (stop - start) + 1 >= ARGS["windowSize"]:

                #sliding window of size win_size
                for i in range(start, stop - ARGS["windowSize"] + 1):

                    dup_seg_sp = synt.to_dup_segments(table_entries_sp[i:i+ARGS["windowSize"]])

                    ind_chrom_one_gene = np.where(np.sum(dup_seg_sp.matrix, axis=0) == 1)[0]

                    for pos in dup_seg_sp.genes_dict:

                        for curr_chrom in dup_seg_sp.genes_dict[pos]:

                            if curr_chrom not in ind_chrom_one_gene:

                                DGENES[spec]["synteny"].update(set(dup_seg_sp.genes_dict[pos]\
                                                                                 [curr_chrom]))

                            else:
                                DGENES[spec]["nosynteny"].update(set(dup_seg_sp.genes_dict[pos]\
                                                                                   [curr_chrom]))
    NO_SYNT = []
    for spec in DGENES:
        DGENES[spec]["nosynteny"] = DGENES[spec]["nosynteny"].difference(DGENES[spec]["synteny"])
        print(spec, len(DGENES[spec]["nosynteny"]))
        NO_SYNT += list(DGENES[spec]["nosynteny"])


    ut.write_updated_orthotable(ALL_ORTHOS, OUTGR, SPECIES, CHROMOSOMES, ARGS["out"],
                                wsize=ARGS['windowSize'], filt_genes=NO_SYNT)
