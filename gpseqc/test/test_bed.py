# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: bed-related method tests.
'''

# DEPENDENCIES =================================================================

import numpy as np
import os
import pandas as pd
import pybedtools as pbt

from gpseqc import bed

# PARAMS =======================================================================

# Bed file, with names and value, bin size 50, bin step 25, chr-size 100
bedstr  = "1\t0\t50\tr1\t1\n"
bedstr += "1\t25\t75\tr2\t2\n"
bedstr += "1\t50\t100\tr3\t4\n"
bedstr += "1\t75\t125\tr4\t5\n"

# Expected windows from bedstr, with {1:100} chr-size
winstr  = "1\t0\t50\n"
winstr += "1\t25\t75\n"
winstr += "1\t50\t100\n"
winstr += "1\t75\t125\n"

# Bins with bin size == bin step = 50, chr-size = 200
nobinstr  = "1\t0\t50\n"
nobinstr += "1\t50\t100\n"
nobinstr += "1\t100\t150\n"
nobinstr += "1\t150\t200\n"

# Bins with bin size = 50, bin step = 25, chr-size = 150
obinstr  = "1\t0\t50\n"
obinstr += "1\t25\t75\n"
obinstr += "1\t50\t100\n"
obinstr += "1\t75\t125\n"
obinstr += "1\t100\t150\n"

# Expected binning of bedstr with nobinstr
nobinnedstr  = "1\t0\t50\trow_1\t1\n"
nobinnedstr += "1\t0\t50\trow_2\t2\n"
nobinnedstr += "1\t50\t100\trow_3\t4\n"
nobinnedstr += "1\t50\t100\trow_4\t5\n"

# Expected binning of bedstr with obinstr
obinnedstr  = "1\t0\t50\trow_1\t1\n"
obinnedstr += "1\t0\t50\trow_2\t2\n"
obinnedstr += "1\t25\t75\trow_3\t1\n"
obinnedstr += "1\t25\t75\trow_4\t2\n"
obinnedstr += "1\t25\t75\trow_5\t4\n"
obinnedstr += "1\t50\t100\trow_6\t2\n"
obinnedstr += "1\t50\t100\trow_7\t4\n"
obinnedstr += "1\t50\t100\trow_8\t5\n"
obinnedstr += "1\t75\t125\trow_9\t4\n"
obinnedstr += "1\t75\t125\trow_10\t5\n"
obinnedstr += "1\t100\t150\trow_11\t5\n"

# Expected combined (sum) binning of bedstr with nobinstr
nobinnedsumstr  = "1\t0\t50\trow_1\t3\n"
nobinnedsumstr += "1\t50\t100\trow_2\t9\n"

# Expected combined (sum) binning of bedstr with obinstr
obinnedsumstr = "1\t0\t50\trow_1\t3\n"
obinnedsumstr += "1\t25\t75\trow_2\t7\n"
obinnedsumstr += "1\t50\t100\trow_3\t11\n"
obinnedsumstr += "1\t75\t125\trow_4\t9\n"
obinnedsumstr += "1\t100\t150\trow_5\t5\n"

# Normalization of bedstr against itself
normstr  = "1\t0\t50\tr1\t1.00\n"
normstr += "1\t25\t75\tr2\t1.00\n"
normstr += "1\t50\t100\tr3\t1.00\n"
normstr += "1\t75\t125\tr4\t1.00\n"

# Statistics for obinnedstr
binstatsdf = pd.DataFrame([
    ['1', 0, 50, 3.0, 1.5, 0.5, 2],
    ['1', 25, 75, 7.0, 7/3., np.std([1, 2, 4]), 3],
    ['1', 50, 100, 11.0, 11/3, np.std([2, 4, 5]), 3],
    ['1', 75, 125, 9.0, 4.5, 0.5, 2],
    ['1', 100, 150, 5.0, 5.0, 0.0, 1]
])
binstatsdf.columns = ["chrom", "start", "end", "sum", "mean", "std", "count"]

# FUNCTIONS ====================================================================

def test_calcStats():
    b = pbt.BedTool(obinnedstr, from_string = True)
    b = bed.calc_stats(b)
    assert b.shape[0] == binstatsdf.shape[0]
    for i in range(b.shape[0]):
        assert b.ix[i, :].tolist() == binstatsdf.ix[i, :].tolist()

def test_mkWindows():
    b = bed.mk_windows({1:100}, 50, 25)
    with open(b.fn, "r+") as IH: content = "".join(IH.readlines())
    assert content == winstr

def test_getChrSize():
    b = pbt.BedTool(bedstr, from_string = True)
    assert bed.get_chr_size(b.fn) == {"1" : 125}

def test_isOverlapping_True():
    assert bed.is_overlapping(pbt.BedTool(obinstr, from_string = True))

def test_isOverlapping_False():
    assert not bed.is_overlapping(pbt.BedTool(nobinstr, from_string = True))

def test_normalize():
    b = pbt.BedTool(bedstr, from_string = True)
    b = bed.normalize(b, b)
    with open(b.fn, "r+") as IH: content = "".join(IH.readlines())
    assert content == normstr

def test_toBins_nonOverlapping():
    a = pbt.BedTool(nobinstr, from_string = True)
    b = pbt.BedTool(bedstr, from_string = True)
    binned = bed.to_bins(a, b)
    with open(binned.fn, "r+") as IH: content = "".join(IH.readlines())
    assert content == nobinnedstr

def test_toBins_overlapping():
    a = pbt.BedTool(obinstr, from_string = True)
    b = pbt.BedTool(bedstr, from_string = True)
    binned = bed.to_bins(a, b)
    with open(binned.fn, "r+") as IH: content = "".join(IH.readlines())
    assert content == obinnedstr

def test_toCombinedBins_nonOverlapping():
    a = pbt.BedTool(nobinstr, from_string = True)
    b = pbt.BedTool(bedstr, from_string = True)
    binned = bed.to_combined_bins(a, b)
    with open(binned.fn, "r+") as IH: content = "".join(IH.readlines())
    assert content == nobinnedsumstr

def test_toCombinedBins_overlapping():
    a = pbt.BedTool(obinstr, from_string = True)
    b = pbt.BedTool(bedstr, from_string = True)
    binned = bed.to_combined_bins(a, b)
    with open(binned.fn, "r+") as IH: content = "".join(IH.readlines())
    assert content == obinnedsumstr

# END ==========================================================================

################################################################################
