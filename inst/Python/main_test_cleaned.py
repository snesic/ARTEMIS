import pandas as pd
import math as math
from init import *
from score import *
from align import *


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

# Impute score matrix, retrieve relevant vars
TSW_scoreMat_cleaned(s1, s1_len, s2, s2_len, g, T, H, TR, TC, traceMat, s, method)

print("H =\n", H)
print("TR =\n", TR)
print("TC =\n", TC)
print("traceMat =\n", traceMat)


