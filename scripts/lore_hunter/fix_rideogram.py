#!/usr/bin/env python

"""
    Fix RIdeogram karyotype figure by adding a legend and title to it.

    Example:

        $ python -m scripts.lore_hunter.fix_rideograms -i fig.svg -o out.svg -l AORe LOre [-c 2]
                                                       [-t '']
"""

import argparse
import os

import matplotlib.pyplot as plt

import seaborn as sns

import svgutils.transform as st


def make_legend(outfilename, title, colors, labels):

    """
    Plots to file a matplotlib figure with only legend and title.

    Args:
        outfilename (str): name for the output figure file
        title (str): title for the figure
        colors (list): ordered list of colors for the legend
        labels (list): ordered list of labels for the legend

    """

    plt.figure(figsize=(8.44, 6))
    func = lambda c: plt.plot([], [], marker='s', color=c, ls="none")[0]
    handles = [func(colors[i]) for i in range(len(colors))]

    plt.axis('off')
    plt.title(title)
    legend = plt.legend(handles, labels, loc="upper right")

    fig = legend.figure
    fig.canvas.draw()
    fig.savefig(outfilename, transparent=True)
    plt.close("all")

def add_legend(input_svg, legend_svg, outfile):
    """
    Create a new svg by putting one svg on top of another.

    Args:
        input_svg (str): name for first svg file
        legend_svg (str): name for second svg file (will be drawn on top of first)
        outfile (list): name for the output figure file

    """
    template = st.fromfile(input_svg)
    second_svg = st.fromfile(legend_svg)
    template.append(second_svg)
    template.save(outfile)

if __name__ == '__main__':

    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required
    PARSER.add_argument('-i', '--input', required=True)

    PARSER.add_argument('-o', '--output', required=True)

    PARSER.add_argument('-l', '--labels', required=True, nargs='+')

    PARSER.add_argument('-c', '--cat', required=False, type=int, default=2)

    PARSER.add_argument('-t', '--title', required=False, default="")


    ARGS = vars(PARSER.parse_args())

    COLORS = sns.color_palette("tab10", ARGS["cat"])

    if ARGS["cat"] == 1:
        COLORS = ["#FF0000"]


    OUT, _ = os.path.splitext(ARGS["output"])

    LEGEND_FILE = OUT + '_legend.svg'

    assert len(ARGS["labels"]) == ARGS["cat"], "Error: different nb of labels and colors in legend."

    make_legend(LEGEND_FILE, ARGS["title"], COLORS, ARGS["labels"])

    add_legend(ARGS["input"], LEGEND_FILE, ARGS["output"])
