"""
    Script to draw and color a genome according to an input classification.
    For instance, it can be used to draw post-duplication chromosomes on modern species.
    The classification provided in input should give a color or class to genes of the species to
    draw. Several input format and plot parameters can be given.

    Example:
        $ python plot_paralogy_map.py -c colored_ancgenes.tsv -g genes.Oryzias.latipes.list.bz2
                                      -o Medaka_paralogymap.svg [-s Medaka] [-f bed] [-maxC 30]
                                      [-minL 30]
"""

# python src/plot_paralogy_map.py -c ../SCORPiOs/ancgenes_inconsistent_trees.tsv
# -g Salmo_salar_genes_chrom_names.bz2 -o test_salmon.svg -s Salmo.salar -sort "names" -f dyogen
# -t "4R sequence/synteny conflicts --save

import sys
import os
import argparse
import pickle
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.colors import is_color_like

import seaborn as sns

from order_chrom import ORDER_CHROM
from palette import PALETTE, REORDER_CHROMS

from scripts.synteny.mygenome import Genome


def read_ancgenes_colors(file_anc_colors, genes, anc=False, species='', out=''):

    """
    Loads predicted colors (PREDUP_CHROM+'_'+POSTDUP_LETTER) for all genes of the considered
    modern species. If `anc` is True load only ancestral genes.


    Args:

        file_anc_colors (str): path to the tab-delimited annotated ancgenes file

        genes (set): set with names of all gene names of the considered species

        anc (bool, optional): whether only ancestral genes should be loaded
                              (if True `genes` can be empty)

        species (str, optional): species name, optional, used to print annotation stats


    Returns:

        dict: gives for each modern gene (key) with prediction its predicted post-duplication
              ancestral chromosome (value).
    """

    assert genes or anc, "empty genes list please check arguments"

    genes_anc = {}
    tot = len(genes)
    pred = 0

    with open(file_anc_colors, 'r') as infile:

        for line in infile:

            line = line.strip().split("\t")

            ancgene, descendants = line[:2]
            anc_chr = line[-1]
            descendants = descendants.split()

            if anc:
                if anc_chr != '?':
                    genes_anc[ancgene] = anc_chr
                    continue

            if anc_chr != '?':

                target_genes = genes.intersection(descendants)

                for gene in target_genes:
                    genes_anc[gene] = anc_chr
                    pred += 1

    if anc:
        pred = len(genes_anc)

    if genes:
        frac = pred/float(tot) * 100
        stat = species+': '+str(pred)+' annotated genes ('+str(round(frac, 2))+'%)\n'

        if out:
            with open(out, 'w') as outfile:
                outfile.write(stat)
        else:
            print(stat)
    return genes_anc


def draw_colors(dgenes, order, genes_colors, species, out, palette=None, min_length=30, max_chr=30,\
                sort_by="size", title='Paralogy Map', save=False):

    """
    Uses matplotlib to draw the genome annotated by post-duplication chromosomes.


    Args:
        dgenes (Genome.genes_list): genome to plot

        order (dict): pre-assigned chromosomes order based on Figures in Nakatani and McLysaght

        genes_colors (dict): for each gene (key) its predicted post-duplication chromosomes (value)

        out (str): path for output figure

        palette (dict): for each post-duplication chromosomes (key) its associated color (value)

        min_length (int, optional): minimum number of genes to plot a chromosome

        max_chr (int, optional): maximum numbre of chromosomes to plot

        sort_by (str, optional): 'size' if chromosomes are to be sorted by size, 'name' by names.

        title (str, optional): title for the generated figure

        save (bool, optional): whether to pickle dump the plotted python dict
    """
    default_palette = sns.color_palette("Set2")
    assert sort_by in ["size", "names"], "Invalid `sort_by` argument, please check"

    #loaded pre-defined chrom order, if specified
    if species in order:
        order = order[species]

    #compute chrom order by size
    elif sort_by == "size":
        order = [str(i) for i in sorted(dgenes.keys(), key=lambda chrom: len(dgenes[chrom]),\
                reverse=True)]

    #compute chrom order by name
    elif sort_by == "names":
        order = sorted([i for i in dgenes.keys() if isinstance(i, int)])
        if not order:
            order = sorted([i for i in dgenes.keys() if i])
        order = [str(i) for i in order]

    i = 0
    chrom_draw = {}
    height = 0.9
    spacing = 0.9
    xranges, colors = [], []
    j = 0

    #load pre-defined color palette
    if not palette:
        palette = {}
        for j, color in enumerate(default_palette):
            palette[str(j)] = color

    #Fill the python dict for plot
    for chromosome in dgenes:
        for gene in dgenes[chromosome]:
            chrom = str(chromosome)
            chrom = chrom.replace("group", "")
            name = gene.names[0]
            col = 'whitesmoke'
            if name in genes_colors:

                col = ''.join(genes_colors[name])
                assert col in palette or is_color_like(col), f'Cannot understand color {col}'

                if col in palette:
                    col = palette[col]

            xranges.append((i, 1))
            colors.append(col)
            i += 1

        chrom_draw[chrom] = (xranges, colors)
        xranges, colors = [], []
        i = 0

    #plot
    _, ax = plt.subplots(1, 1)
    plt.title(species +" "+title)
    yticks = []
    yticklabels = []
    ymin, nb = 0, 0
    to_save = {}
    for chrom in order:
        xranges, colors = chrom_draw[chrom]
        to_save[chrom] = colors
        if len(colors) > min_length and nb < max_chr:
            ymin += height + spacing
            yrange = (ymin, height)
            coll = BrokenBarHCollection(xranges, yrange, facecolors=colors)
            ax.add_collection(coll)
            center = yrange[0] + yrange[1]/2.
            yticks.append(center)
            yticklabels.append(chrom)
            nb += 1

    ax.axis('tight')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    plt.ylabel(species+ ' chromosomes')
    plt.xlabel("Genomic position (in genes unit)")
    sns.despine()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close('all')

    #dump dict to file
    if save:

        #for recombination pipeline
        # outpkl, _ = os.path.splitext(out)
        # with open(outpkl+'.pkl', "wb") as outf:

        #for paralogy_map_pipeline
        with open(species+"_"+title+'.pkl', "wb") as outf:

            pickle.dump(to_save, outf)




def draw_anc(dgenes, order, species, out, palette, min_length=30, max_chr=30,\
             sort_by="size", title='Paralogy Map', save=False):

    """
    Uses matplotlib to draw the genome annotated by post-duplication chromosomes.


    Args:
        dgenes (Genome.genes_list): genome to plot

        order (dict): pre-assigned chromosomes order based on Figures in Nakatani and McLysaght

        genes_colors (dict): for each gene (key) its predicted post-duplication chromosomes (value)

        out (str): path for output figure

        palette (dict): for each post-duplication chromosomes (key) its associated color (value)

        min_length (int, optional): minimum number of genes to plot a chromosome

        max_chr (int, optional): maximum numbre of chromosomes to plot

        sort_by (str, optional): 'size' if chromosomes are to be sorted by size, 'name' by names.

        title (str, optional): title for the generated figure

        save (bool, optional): whether to pickle dump the plotted python dict
    """

    default_palette = sns.color_palette("Set2")
    assert sort_by in ["size", "names"], "Invalid `sort_by` argument, please check"

    #loaded pre-defined chrom order, if specified
    if species in order:
        order = order[species]

    #compute chrom order by size
    elif sort_by == "size":
        order = [str(i) for i in sorted(dgenes.keys(), key=lambda chrom: len(dgenes[chrom]),\
                reverse=True)]

    #compute chrom order by name
    elif sort_by == "names":
        order = sorted([i for i in dgenes.keys() if isinstance(i, int)])
        order = [str(i) for i in order]

    i = 0
    chrom_draw = {}
    height = 0.9
    spacing = 0.9
    xranges, colors = [], []
    j = 0

    #load pre-defined color palette
    if not palette:
        for j, color in enumerate(default_palette):
            palette[str(j)] = color

    order = []
    #Fill the python dict for plot
    for chromosome in range(1, 14):
        for letter in ["a", "b"]:
            chrom = str(chromosome) + letter
            order.append(chrom)
            for _ in dgenes[chrom]:

                col = chrom
                assert col in palette or is_color_like(col), f'Cannot understand color {col}'

                if col in palette:
                    col = palette[col]

                xranges.append((i, 1))
                colors.append(col)
                i += 1

            chrom_draw[chrom] = (xranges, colors)
            xranges, colors = [], []
            i = 0

    #plot
    _, ax = plt.subplots(1, 1)
    plt.title(species +" "+title)
    yticks = []
    yticklabels = []
    ymin, nb = 0, 0
    to_save = {}
    for chrom in order[::-1]:
        xranges, colors = chrom_draw[chrom]
        to_save[chrom] = colors
        if len(colors) > min_length and nb < max_chr:
            ymin += height + spacing
            yrange = (ymin, height)
            coll = BrokenBarHCollection(xranges, yrange, facecolors=colors)
            ax.add_collection(coll)
            center = yrange[0] + yrange[1]/2.
            yticks.append(center)
            yticklabels.append(chrom)
            nb += 1

    ax.axis('tight')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    plt.ylabel('ancestral chromosomes')
    plt.xlabel("Number of genes")
    sns.despine()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close('all')

    #dump dict to file
    if save:
        with open(species+"_"+title+'.pkl', "wb") as outf:
            pickle.dump(to_save, outf)

def load_palette_from_file(input_file):

    """
    Loads a color palette from file.

    Args:
        input_file (str): Input file with chromosome name and associated color

    Returns:
        dict: for each chromosome (key) its associated color for plot (value)
    """

    palette = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            chrom, col = line.strip().split("\t")
            palette[chrom] = col
    return palette


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    ## Required ##
    PARSER.add_argument('-c', '--color', type=str, help='Input file with classification',
                        required=True)

    PARSER.add_argument('-g', '--genes', type=str, help='Genes file of the target species',
                        required=False, default='')


    ## Optional ##
    PARSER.add_argument('-o', '--output_file', type=str, help='output file', required=False,
                        default='out.svg')

    PARSER.add_argument('-s', '--species_name', type=str, help="Name of the species for plot title",
                        required=False, default='')

    PARSER.add_argument('-f', '--genesformat', type=str, required=False, help="Format of the genes"
                        "file either bed or dyogen", default='bed')

    PARSER.add_argument('-maxC', '--max_chr', type=int, required=False, help="Max number of\
                        chromosomes to plot", default=30)

    PARSER.add_argument('-minL', '--min_length', type=int, required=False, help="Minimum length for\
                        a chromosome or scaffold to be plotted", default=30)

    PARSER.add_argument('-sort', '--sort_by', type=str, required=False, help="Order for chromosomes\
                        can be 'size' or 'names'", default="size")

    PARSER.add_argument('-t', '--title', type=str, required=False, help="Plot title, species name\
                        will be appended, if provided", default="Paralogy Map")

    PARSER.add_argument('-pf', '--palette_from_file', type=str, help="Redefine color palette wth a"
                        " tab-delimited file, giving a class-color correspondence",
                        required=False, default='')

    PARSER.add_argument('--save', action='store_true', help="If specified, pickles the dictionnary\
                        used for plot")

    PARSER.add_argument('--singlesp', action='store_true', help="Specifies that input is not in \
                        ancGenes format but contains only genes of the target species.")

    PARSER.add_argument('-gid', '--gene_index', type=int, required=False, default=0, help="Use in\
                        conjonction with --singlesp, specify column with gene names.\
                        Classes are always expected to be the last column")

    PARSER.add_argument('-os', '--stats_out', type=str, required=False, default='', help="File to\
                        dump stats")

    PARSER.add_argument('--dontdraw', action='store_true')

    PARSER.add_argument('--default_palette', action='store_true')

    ARGS = vars(PARSER.parse_args())

    if not ARGS["palette_from_file"] and not ARGS["default_palette"]:

        #use color of "a" genes to serve as colors for pre-duplication chr
        KEYS = set(PALETTE.keys())
        for key in KEYS:
            if key[-1] == 'a' and key != '9a':
                PALETTE[key[:-1]] = PALETTE[key]
            elif key == '9b':
                PALETTE[key[:-1]] = PALETTE[key]
    elif ARGS["palette_from_file"]:

        PALETTE = load_palette_from_file(ARGS["palette_from_file"])

    else:
        PALETTE = None

    GENES = {}

    if ARGS["genes"]:

        GENOME = Genome(ARGS["genes"], ARGS["genesformat"])
        GENES = {g.names[0] for g in GENOME}

        GENES_COL = read_ancgenes_colors(ARGS["color"], GENES, anc=ARGS['singlesp'],
                                         species=ARGS["species_name"], out=ARGS["stats_out"])

        if not ARGS["dontdraw"]:
            draw_colors(GENOME.genes_list, ORDER_CHROM, GENES_COL, ARGS["species_name"],
                        ARGS["output_file"], PALETTE, min_length=ARGS["min_length"],
                        max_chr=ARGS["max_chr"], title=ARGS["title"], sort_by=ARGS["sort_by"],
                        save=ARGS['save'])

        else:
            sys.stderr.write("Touching empty figure file as drawing is disabled with --dontdraw\n")
            Path(ARGS["output_file"]).touch(exist_ok=True)

    else:

        GENES_COL = read_ancgenes_colors(ARGS["color"], GENES, anc=ARGS['singlesp'],
                                         species=ARGS["species_name"], out=ARGS["stats_out"])


        PALETTE_NEW = {}
        for val in GENES_COL.values():
            if val != '?':
                new_val = str(REORDER_CHROMS[int(val[:-1])]) + val[-1]
                GENES[new_val] = {i for i in GENES_COL if GENES_COL[i] == val}
                PALETTE_NEW[new_val] = PALETTE[val]
        draw_anc(GENES, ORDER_CHROM, ARGS["species_name"],
                 ARGS["output_file"], PALETTE_NEW, min_length=ARGS["min_length"],
                 max_chr=ARGS["max_chr"], title=ARGS["title"], sort_by=ARGS["sort_by"],
                 save=ARGS['save'])

        # PALETTE_NEW = {}
        # seen = []
        # for val in GENES_COL.values():
        #     if val != '?':
        #         val = int(val[:-1])
        #         if val not in seen:
        #             print(val)
        #             seen.append(val)

        #             # new_val = str(REORDER_CHROMS[int(val[:-1])]) + val[-1]
        #             GENES[str(REORDER_CHROMS[val])+'b'] = {'_'.join(i.split('_')[:-1]) for i in GENES_COL if GENES_COL[i][:-1] and int(GENES_COL[i][:-1])==val}
        #             GENES[str(REORDER_CHROMS[val])+'a'] = GENES[str(REORDER_CHROMS[val])+'b']
        #             PALETTE_NEW[str(REORDER_CHROMS[val])+'a'] = PALETTE[str(val)+'a']
        #             PALETTE_NEW[str(REORDER_CHROMS[val])+'b'] = PALETTE[str(val)+'b']
        # draw_anc(GENES, ORDER_CHROM, ARGS["species_name"],
        #          ARGS["output_file"], PALETTE_NEW, min_length=ARGS["min_length"],
        #          max_chr=ARGS["max_chr"], title=ARGS["title"], sort_by=ARGS["sort_by"],
        #          save=ARGS['save'])
    #TODO: option to draw pre and post duplication chromosomes -
    #-> generalized funtion instead of ugly hack above
