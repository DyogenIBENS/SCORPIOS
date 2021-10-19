#!/usr/bin/env python

"""
    This script loads 2 scores distributions and finds the optimal discriminative threshold to
    separate distributions based on the F1-score, assuming true positives to recover are in the
    distribution of higher scores.

    Inputs are python lists pickled in files, output is written to file with the :code:`--support`
    prefix, to call the script missed_orthologies.py in snakemake with the :code:`--support` arg.

    Example::

        $ python -m scripts.synteny.f1_score_optimization -i1 scores_1.pkl -i2 scores_2.pkl
        [-out out]
"""


import argparse
import pickle

import numpy as np

def load_scores(input1, input2):

    """
    Unpickles the lists of scores.

    Args:
        input1 (str) : paths to the pickled object 1
        input2 (str) : paths to the pickled object 2

    Returns:

        tuple: a tuple containing:

            scores1, scores2: the unpickled lists
    """

    with open(input1, 'rb') as infile1, open(input2, 'rb') as infile2:
        scores1 = pickle.load(infile1)
        scores2 = pickle.load(infile2)

    if np.median(scores1) < np.median(scores2):
        scores1, scores2 = scores2, scores1

    return scores1, scores2


def compute_f1(scores1, scores2, threshold):

    """
    Computes the F1-score for a given threshold.

    Args:
        scores1 (list): list of scores 1
        scores2 (list): list of scores 2
        threshold (float): threshold value

    Returns:
        float: F1-score
    """

    true_pos = len([x for x in scores1 if x >= threshold])
    false_neg = len(scores1) - true_pos
    false_pos = len([x for x in scores2 if x >= threshold])

    recall = true_pos / float(true_pos + false_neg) if (true_pos + false_neg) else 0
    precision = true_pos / float(true_pos + false_pos) if (true_pos + false_pos) else 0

    f1_score = 2*precision*recall / (precision + recall) if (precision + recall) else 0

    return f1_score


def get_discriminant_threshold(input1, input2, test_range=[j for j in range(30)]):

    """
    Finds the most discriminative threshold between the two distributions based on F1-score.

    Args:
        input1, input2 (str): paths to the pickled objects
        test_range (list, optional): list of thresholds to test

    Returns:
        int: optimized threshold based on F1-score
    """

    scores1, scores2 = load_scores(input1, input2)

    if not scores1 or not scores2:
        return 2

    best_threshold = 0
    max_f1 = 0
    for i in test_range:
        f1_score = compute_f1(scores1, scores2, i)
        f1_score = round(f1_score, 2) #take the most conservative threshold for improvement < 0.01
        if f1_score >= max_f1:
            best_threshold = i
            max_f1 = f1_score
    return best_threshold


if __name__ == '__main__':

    #Arguments

    PARSER = argparse.ArgumentParser(description=__doc__,\
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    #Required

    PARSER.add_argument('-i1', '--input1', help='First score distribution', required=True)

    PARSER.add_argument('-i2', '--input2', help='Second score distribution', required=True)


    #Optional

    PARSER.add_argument('-out', '--output', help='Output file.', required=False, default="out")

    ARGUMENTS = vars(PARSER.parse_args())

    THRESHOLD = get_discriminant_threshold(ARGUMENTS['input1'], ARGUMENTS["input2"])

    with open(ARGUMENTS["output"], 'w') as OUTFILE:

        if THRESHOLD:
            OUTFILE.write('--support '+str(round(THRESHOLD, 1)))
