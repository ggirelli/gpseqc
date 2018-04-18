# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: bed-related methods.
'''

# DEPENDENCIES =================================================================

import numpy as np
import os
import pandas as pd
import pybedtools as pbt

from gpseq_ce import bed

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

# Bins with bin size == bin step = 50, chr-size = 100
binstr  = "1\t0\t50\n"
binstr += "1\t50\t100\n"

# Expected binning of bedstr with binstr
binnedstr  = "1\t0\t50\trow_1\t1\n"
binnedstr += "1\t0\t50\trow_2\t2\n"
binnedstr += "1\t50\t100\trow_3\t4\n"
binnedstr += "1\t50\t100\trow_4\t5\n"

# Expected combined (sum) binning of bedstr with binstr
binnedsumstr  = "1\t0\t50\trow_1\t3\n"
binnedsumstr += "1\t50\t100\trow_2\t9\n"

# FUNCTIONS ====================================================================

def test_calcStats():
    pass

def test_mkWindows():
    b = bed.mk_windows({1:100}, 50, 25)
    with open(b.fn, "r+") as IH:
        content = "".join(IH.readlines())
    assert content == winstr

def test_getChrSize():
    b = pbt.BedTool(bedstr, from_string = True)
    assert bed.get_chr_size(b.fn) == {"1" : 125}

def test_normalize():
    pass

def test_toBins():
    a = pbt.BedTool(binstr, from_string = True)
    b = pbt.BedTool(bedstr, from_string = True)
    binned = bed.to_bins(a, b)
    with open(binned.fn, "r+") as IH:
        content = "".join(IH.readlines())
    assert content == binnedstr

def test_toCombinedBins():
    a = pbt.BedTool(binstr, from_string = True)
    b = pbt.BedTool(bedstr, from_string = True)
    binned = bed.to_combined_bins(a, b)
    with open(binned.fn, "r+") as IH:
        content = "".join(IH.readlines())
    assert content == binnedsumstr

# END ==========================================================================

################################################################################
