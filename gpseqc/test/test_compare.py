# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: ranking comparison-related method tests.
'''

# DEPENDENCIES =================================================================

import numpy as np
import pandas as pd
from scipy.stats import norm

from gpseqc.compare import *

# PARAMS =======================================================================

r1 = pd.DataFrame.from_items([
    ("chrom",    ["A", "B", "C", "D"]),
    ("start",    ["0", "0", "0", "0"]),
    ("end",      ["0", "0", "0", "0"]),
    ("prob_2p",  [ 1,   3,   5,   20]),
    ("prob_g",   [ 6,   20,  5,   7 ])
])
rt1 = RankTable(df = r1)

r2 = pd.DataFrame.from_items([
    ("chrom",    ["A", "C", "D"]),
    ("start",    ["0", "0", "0"]),
    ("end",      ["0", "0", "0"]),
    ("prob_2p",  [ 1,   5,   20]),
    ("prob_g",   [ 6,   5,   7 ])
])
rt2 = RankTable(df = r2)

r3 = pd.DataFrame.from_items([
    ("chrom",    ["A", "C", "D", "N"]),
    ("start",    ["0", "0", "0", "0"]),
    ("end",      ["0", "0", "0", "0"]),
    ("prob_2p",  [ 1,   5,   20, np.nan]),
    ("prob_g",   [ 6,   5,   7,  5 ])
])
rt3 = RankTable(df = r3)

rs = [("A", "0", "0"), ("B", "0", "0"), ("C", "0", "0"), ("D", "0", "0")]

# FUNCTIONS ====================================================================

def test_RankTable():
    assert all(r1.iloc[:, :4] == rt1[0]._df) # 
    assert 2 == len(rt1)                     # Number of metrics
    assert 4 == len(rt1[0])                  # Number of regions
    assert ["prob_2p", "prob_g"] == [m for m in rt1.available_metrics()]
    assert rs == rt1._all_regions()          # Regions
    assert [rs[0], rs[2], rs[3]] == (rt1 & rt2)._all_regions()
    assert len(rt1.shuffle()._all_regions()) == len(rs)
    assert all([r in rt1.shuffle()._all_regions() for r in rs])

def test_MetricTable():
    mt1 = rt1[0]
    mt2 = rt1[1]
    mt3 = rt2[0]
    assert 4 == len(mt1)
    assert "prob_2p" == mt1._metric
    assert ("A", "0", "0", 1) == tuple(mt1[0].tolist())
    assert ("C", "0", "0", 5) == tuple(mt2[0].tolist())
    assert (1, 3, 5, 20) == tuple(mt1.mcol.tolist())
    assert (5, 6, 7, 20) == tuple(mt2.mcol.tolist())
    assert [rs[0], rs[2], rs[3]] == (mt1 & mt3)._all_regions()
    assert len(mt1.shuffle()._all_regions()) == len(rs)

def test_MetricTable_KendallTau():
    mt1 = rt1[0]
    mt2 = rt1[1]
    mt3 = rt3[0]
    assert 0 == mt1.calc_KendallTau(mt1)
    assert 0.5 == mt1.calc_KendallTau(mt2)
    assert 0 == mt1.calc_KendallTau(mt3)
    assert 1/3. == mt2.calc_KendallTau(mt3)

def test_MetricTable_KendallTau_weighted():
    mt1 = rt1[0]
    mt2 = rt1[1]
    mt3 = rt3[0]
    assert 0 == mt1.calc_KendallTau_weighted(mt1)
    assert 0.5101 == np.round(dKTw_iter(0, 1, rt1, rt1), 4)
    assert 0.5101 == np.round(mt1.calc_KendallTau_weighted(mt2), 4)
    assert 0 == mt1.calc_KendallTau(mt3)
    assert 0.1776 == np.round(mt2.calc_KendallTau_weighted(mt3), 4)

def test_compare2randDistr():
    np.random.seed(654546)
    g = np.random.randn(5000) * 5
    mu, sigma = norm.fit(g)
    Zpval = norm.cdf(0, mu, sigma)
    if .5 < Zpval: Zpval = 1 - Zpval
    Zpval *= 2
    assert 0.9822733103866876 == Zpval

def test_mk2DdistanceMatrix():
    distance_matrix = mk2DdistanceMatrix(10, 10)
    for i in range(distance_matrix.shape[0]):
        for j in range(distance_matrix.shape[0]):
            assert np.absolute(j - i) == distance_matrix[i, j]

def test_EMD():
    a_weights = r1.iloc[:, 3]
    a_weights /= a_weights.sum()
    b_weights = r1.iloc[:, 4]
    b_weights /= b_weights.sum()

    distance_matrix = mk2DdistanceMatrix(len(a_weights), len(b_weights))

    a_isorted = np.argsort(a_weights)
    b_asorted = b_weights[a_isorted]
    a_asorted = a_weights[a_isorted]
    d1  = calc_EarthMoversDistance(a_asorted, b_asorted, distance_matrix)
    assert 0.39171 == np.round(d1, 5)

    b_isorted = np.argsort(b_weights)
    a_bsorted = a_weights[b_isorted]
    b_bsorted = b_weights[b_isorted]
    d2 = calc_EarthMoversDistance(a_bsorted, b_bsorted, distance_matrix)
    assert 0.18209 == np.round(d2, 5)

    a_extreme = np.zeros(a_weights.shape)
    a_extreme[0] = a_weights.sum()
    b_extreme = np.zeros(b_weights.shape)
    b_extreme[-1] = b_weights.sum()
    assert distance_matrix.max() == np.round(pyemd.emd(a_extreme, b_extreme,
        distance_matrix, extra_mass_penalty = -1.0))

    assert 0.28690 == np.round(dEMD_iter(0, 1, rt1, rt1), 5)

    t1 = rt1 & rt3
    t3 = rt3 & rt1
    assert 0 == np.round(dEMD_iter(0, 0, t1, t3), 5)
    assert 0.28526 == np.round(dEMD_iter(1, 0, t1, t3), 5)

# END ==========================================================================

################################################################################
