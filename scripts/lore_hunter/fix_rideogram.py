import argparse
import os

import svgutils.transform as st

import matplotlib.pyplot as plt

import seaborn as sns


def make_legend(outfilename, title, colors):
    """
    """
    plt.figure(figsize=(8.44, 6))
    f = lambda m,c: plt.plot([],[],marker=m, color=c, ls="none")[0]
    handles = [f("s", colors[i]) for i in range(len(colors))]
    labels = ["AORe", "LORe"]
    plt.axis('off')
    plt.title(title)
    legend = plt.legend(handles, labels, loc="upper right")

    fig = legend.figure
    fig.canvas.draw()
    fig.savefig(outfilename, transparent=True)
    plt.close("all")

def add_legend(input_svg, legend_svg, outfile):
    """
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

    PARSER.add_argument('-c', '--cat', required=False, type=int, default=2)

    PARSER.add_argument('-t', '--title', required=False, default="")

    ARGS = vars(PARSER.parse_args())

    COLORS = sns.color_palette("tab10", ARGS["cat"])

    OUT, _ = os.path.splitext(ARGS["output"])

    legend_file = OUT + '_legend.svg'

    make_legend(legend_file, ARGS["title"], COLORS)

    add_legend(ARGS["input"], legend_file, ARGS["output"])