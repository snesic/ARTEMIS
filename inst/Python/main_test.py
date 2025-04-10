import numpy as np
import pandas as pd
import math as math
from re import search
from re import findall
from init import *
from score import *
from align import *

pd.options.display.max_columns = None


def find_gaps(pat, seq):
    gaps_init = search(pat, seq)
    if gaps_init is not None:
        gaps = len(findall("__", gaps_init[0]))
    else:
        gaps = 0

    return gaps


# time gap example:
s1 = [
    ["0", "a"],
    ["0", "c"],
    ["0", "b"],
    ["1", "b"],
    ["1", "b"],
]
s2 = [
    ["1", "b"],
    ["1", "b"],
    ["10", "a"],
    ["0", "c"],
    ["0", "b"],
    ["1", "b"],
]

# simple example, double alignment
s1 = [
    ["0", "a"],
    ["0", "b"],
    ["0", "c"],
]
s2 = [
    ["10", "a"],
    ["1", "b"],
    ["0", "c"],
    ["10", "a"],
    ["0", "b"],
    ["0", "c"],
]

# Time issue
# seq = c("0.n;14.n;14.n;14.n;21.n; 14.n;14.n;14.n;14.n;14.n;14.n;14.n;14.n;13.n;"

s1 = [
    ["1", "a"],
    ["1", "a"],
]
s2 = [
    ["0", "a"],
    ["10", "a"],
    ["11", "a"],
]


s = {
    "a": [1.0, -1.1, -1.1],
    "b": [-1.1, 1.0, -1.1],
    "c": [-1.1, -1.1, 1.0],
}
s = pd.DataFrame(s, index=["a", "b", "c"])

regName = "ee"
g = 0.4
T = 0.5
method = "PropDiff"
mem = -1
verbose = True

s1_len = len(s1)
s2_len = len(s2)

# Initialise the 3 score matrices and the traceback matrix
H = init_Hmat(s1_len, s2_len, g)
TR = init_TRmat(s1, s1_len, s2, s2_len)
TC = init_TCmat(s1, s1_len, s2, s2_len)
traceMat = init_traceMat(s1_len, s2_len)

# Track if secondary alignments have been collected
secondary = 0

# Setup pattern for detecting sequence lengths, by number of "."s (Aligned drugs)
pat = "\."

# Setup pattern for detecting sequence gaps, by number of "__"s (Aligned gaps)
pat_gap = "(__;)+[0-9]+"
# pat_gap = "__"
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

from copy import deepcopy

_s1 = deepcopy(s1)
_s1_len = deepcopy(s1_len)
_s2 = deepcopy(s2)
_s2_len = deepcopy(s2_len)
_g = deepcopy(g)
_T = deepcopy(T)
_H = deepcopy(H)
_TR = deepcopy(TR)
_TC = deepcopy(TC)
_traceMat = deepcopy(traceMat)
_s = deepcopy(s)
_method = deepcopy(method)


TSW_scoreMat_cleaned(_s1,_s1_len,_s2,_s2_len,_g,_T,_H,_TR,_TC,_traceMat, _s, _method)

assert np.array_equal(H, _H), "H matrices are not equal"
print(H.shape, _H.shape)
assert np.array_equal(TR, _TR), "TR matrices are not equal"
print(TR.shape, _TR.shape)
assert np.array_equal(TC, _TC), "TC matrices are not equal"
print(TC.shape, _TC.shape)
assert np.array_equal(traceMat, _traceMat), "traceMat matrices are not equal"
print(traceMat.shape, _traceMat.shape)

print("H =\n", H)
print("TR =\n", TR)
print("TC =\n", TC)
print("traceMat =\n", traceMat)


# Find best scoring cell
finalScore, finalIndex, mem_index, mem_score = find_best_score(
    H, s1_len, s2_len, mem, verbose
)


s1_aligned_t, s2_aligned_t, totAligned_t = align_TSW(
    traceMat, s1, s2, s1_len, s2_len, finalIndex
)

s1_aligned_t
s2_aligned_t


s_f_len = max(len(findall(pat, s2_aligned_t)), len(findall(pat, s1_aligned_t)))

s1_gaps = find_gaps(pat_gap, s1_aligned_t)
s1_end_gaps = find_gaps(pat_end_gap, s1_aligned_t)
s2_gaps = find_gaps(pat_gap, s2_aligned_t)
s2_end_gaps = find_gaps(pat_end_gap, s2_aligned_t)

s1_start = finalIndex[1] - s_f_len
s1_end = finalIndex[1]
s2_start = finalIndex[0] - s_f_len - s1_end_gaps
s2_end = finalIndex[0] - s1_end_gaps

s2_start + 1
s1_start + 1
