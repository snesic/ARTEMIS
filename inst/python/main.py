import numpy as np
import pandas as pd
import math as math
from re import search
from re import findall
from numpy import append
from init import init_Hmat, init_TRmat, init_TCmat, init_traceMat
from score import TSW_scoreMat, find_best_score
from align import align_TSW

pd.options.display.max_columns = None


def find_gaps(pat, seq):
    gaps_init = search(pat, seq)
    if gaps_init is not None:
        gaps = len(findall("__", gaps_init[0]))
    else:
        gaps = 0

    return gaps


def temporal_alignment(
    s1, regName, s2, g, T, s, verbose, mem=-1, removeOverlap=0, method="PropDiff"
):
    s1_len = len(s1)
    s2_len = len(s2)

    # Initialise the 3 score matrices and the traceback matrix
    H = init_Hmat(s1_len, s2_len)
    TR = init_TRmat(s1, s1_len, s2, s2_len)
    TC = init_TCmat(s1, s1_len, s2, s2_len)
    traceMat = init_traceMat(s1_len, s2_len)

    # Track if secondary alignments have been collected
    secondary = 0

    # Setup pattern for detecting sequence lengths, by number of "."s (Aligned drugs)
    pat = "\."

    # Setup pattern for detecting sequence gaps, by number of "__"s (Aligned gaps)
    pat_gap = "(__;)+[0-9]+"
    pat_end_gap = "(__;)+__$|__$"
    pat_search = "__"

    # Init return Dat
    returnDat = [
        regName,
        str(s1).strip("[]"),
        str(s2).strip("[]"),
        "",
        "",
        "",
        "",
        "",
        "",
        "",
    ]
    returnDat = np.array(returnDat, dtype=object)

    # Impute score matrix, retrieve relevant vars
    TSW_scoreMat(s1, s1_len, s2, s2_len, g, T, H, TR, TC, traceMat, s, method)

    # Find best scoring cell
    finalScore, finalIndex, mem_index, mem_score = find_best_score(
        H, s1_len, s2_len, mem, verbose
    )

    for i in range(0, len(mem_index)):
        s1_aligned_t, s2_aligned_t, totAligned_t = align_TSW(
            traceMat, s1, s2, s1_len, s2_len, mem_index[i]
        )

        s_f_len = max(len(findall(pat, s2_aligned_t)), len(findall(pat, s1_aligned_t)))

        s1_gaps = find_gaps(pat_gap, s1_aligned_t)
        s1_end_gaps = find_gaps(pat_end_gap, s1_aligned_t)
        s2_gaps = find_gaps(pat_gap, s2_aligned_t)
        s2_end_gaps = find_gaps(pat_end_gap, s2_aligned_t)

        s1_start = mem_index[i][1] - s_f_len
        s1_end = mem_index[i][1]
        s2_start = mem_index[i][0] - s_f_len + s2_gaps + s2_end_gaps
        s2_end = mem_index[i][0] - s1_end_gaps

        if (s1_start + 1) > 1:
            totAligned_t = totAligned_t + (s1_end - (s1_start + 1))
            s_f_len = s_f_len + (s1_end - (s1_start + 1))

        returnDat = append(
            returnDat,
            [
                regName,
                s1_aligned_t,
                s2_aligned_t,
                mem_score[i],
                s1_start + 1,
                s1_end,
                s2_start + 1,
                s2_end,
                s_f_len,
                totAligned_t,
            ],
            axis=0,
        )

    # Reshape return array to account for secondary alignments
    returnDat = returnDat.reshape(len(mem_index) + 1, 10)
    returnDat = pd.DataFrame(returnDat)

    return returnDat
