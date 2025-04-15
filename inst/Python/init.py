import numpy as np
import pandas as pd
import math as math
import re


def init_Hmat(s1_len, s2_len):
    H = np.zeros((s2_len + 1, s1_len + 1), float)
    return H


def init_TRmat(s1, s1_len, s2, s2_len):
    TR = np.zeros((s2_len + 1, s1_len + 1), float)
    return TR


def init_TCmat(s1, s1_len, s2, s2_len):
    TC = np.zeros((s2_len + 1, s1_len + 1), float)
    return TC


def init_traceMat(s1_len, s2_len):
    traceMat = np.zeros((s2_len + 1, s1_len + 1), float)
    return traceMat
