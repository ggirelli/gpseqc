# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: ranking comparison-related method tests.
'''

# DEPENDENCIES =================================================================

import numpy as np
import pandas as pd

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
    assert 0 == mt1.dKT(mt1)
    assert 0.5 == mt1.dKT(mt2)

def test_MetricTable_KendallTau_weighted():
    mt1 = rt1[0]
    mt2 = rt1[1]
    assert 0 == mt1.dKTw(mt1)
    assert 0.5101 == np.round(mt1.dKTw(mt2), 4)

# END ==========================================================================

################################################################################
