# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: utsite-related method tests.
'''

# DEPENDENCIES =================================================================

import numpy as np
import pandas as pd
import pybedtools as pbt

from gpseq_ce import cutsite_domain as csd

# PARAMS =======================================================================

bedstr1  ='chr1\t500\t506\tr1\t3\n'
bedstr1 +='chr1\t600\t606\tr2\t2\n'
bedstr1 +='chr1\t1000\t1006\tr3\t4\n'
bedstr1 +='chr1\t2050\t2056\tr4\t4\n'
bedstr1 +='chr1\t2100\t2106\tr5\t8\n'
bedstr1 +='chr2\t0\t6\tr7\t5\n'
bedstr1 +='chr2\t20\t26\tr8\t2\n'
bedstr1 +='chr2\t10000\t10006\tr9\t2\n'
bedstr1 +='chr2\t12000\t12006\tr11\t4\n'
bed1 = pbt.BedTool(bedstr1, from_string = True)

bedstr2  = 'chr1\t500\t506\tr1\t1\n'
bedstr2 += 'chr1\t600\t606\tr2\t2\n'
bedstr2 += 'chr1\t1000\t1006\tr3\t5\n'
bedstr2 += 'chr1\t2000\t2006\tr4\t2\n'
bedstr2 += 'chr1\t2100\t2106\tr5\t6\n'
bedstr2 += 'chr2\t0\t6\tr7\t4\n'
bedstr2 += 'chr2\t10\t16\tr8\t1\n'
bedstr2 += 'chr2\t10000\t10006\tr9\t1\n'
bedstr2 += 'chr2\t12000\t12006\tr11\t1\n'
bed2 = pbt.BedTool(bedstr2, from_string = True)

csstr  = 'chr1\t500\t506\tr1\n'
csstr += 'chr1\t600\t606\tr2\n'
csstr += 'chr1\t1000\t1006\tr3\n'
csstr += 'chr1\t2000\t2006\tr4\n'
csstr += 'chr1\t2050\t2056\tr5\n'
csstr += 'chr1\t2100\t2106\tr6\n'
csstr += 'chr2\t0\t6\tr8\n'
csstr += 'chr2\t10\t16\tr9\n'
csstr += 'chr2\t20\t26\tr10\n'
csstr += 'chr2\t10000\t10006\tr11\n'
csstr += 'chr2\t12000\t12006\tr13\n'
csstr += 'chr2\t12100\t12106\tr14\n'
csbed = pbt.BedTool(csstr, from_string = True)

gstr  = 'chr1\t0\t1000\n'
gstr += 'chr1\t1000\t2000\n'
gstr += 'chr1\t2000\t3000\n'
gstr += 'chr2\t0\t1000\n'
gstr += 'chr2\t1000\t2000\n'
gstr += 'chr2\t2000\t3000\n'
gstr += 'chr2\t3000\t4000\n'
gstr += 'chr2\t4000\t5000\n'
gstr += 'chr2\t5000\t6000\n'
gstr += 'chr2\t6000\t7000\n'
gstr += 'chr2\t7000\t8000\n'
gstr += 'chr2\t8000\t9000\n'
gstr += 'chr2\t9000\t10000\n'
gstr += 'chr2\t10000\t11000\n'
gstr += 'chr2\t11000\t12000\n'
gstr += 'chr2\t12000\t13000\n'
gbed = pbt.BedTool(gstr, from_string = True)

ngunionstr  = 'chr1\t0\t1000\n'
ngunionstr += 'chr1\t1000\t2000\n'
ngunionstr += 'chr1\t2000\t3000\n'
ngunionstr += 'chr2\t0\t1000\n'
ngunionstr += 'chr2\t10000\t11000\n'
ngunionstr += 'chr2\t12000\t13000\n'

ngintersecstr  = 'chr1\t500\t506\n'
ngintersecstr += 'chr1\t600\t606\n'
ngintersecstr += 'chr1\t1000\t1006\n'
ngintersecstr += 'chr1\t2050\t2056\n'
ngintersecstr += 'chr1\t2100\t2106\n'
ngintersecstr += 'chr2\t0\t6\n'
ngintersecstr += 'chr2\t20\t26\n'
ngintersecstr += 'chr2\t10000\t10006\n'
ngintersecstr += 'chr2\t12000\t12006\n'

gintersecstr  = 'chr1\t0\t1000\n'
gintersecstr += 'chr1\t1000\t2000\n'
gintersecstr += 'chr1\t2000\t3000\n'
gintersecstr += 'chr2\t0\t1000\n'
gintersecstr += 'chr2\t10000\t11000\n'
gintersecstr += 'chr2\t12000\t13000\n'

appUniverseNGstr  = 'chr1\t500\t506\trow_1\t3\n'
appUniverseNGstr += 'chr1\t600\t606\trow_2\t2\n'
appUniverseNGstr += 'chr1\t1000\t1006\trow_3\t4\n'
appUniverseNGstr += 'chr1\t2050\t2056\trow_4\t4\n'
appUniverseNGstr += 'chr1\t2100\t2106\trow_5\t8\n'
appUniverseNGstr += 'chr2\t0\t6\trow_6\t5\n'
appUniverseNGstr += 'chr2\t20\t26\trow_7\t2\n'
appUniverseNGstr += 'chr2\t10000\t10006\trow_8\t2\n'
appUniverseNGstr += 'chr2\t12000\t12006\trow_9\t4\n'

appUniverseGstr  = 'chr1\t0\t1000\trow_1\t5\n'
appUniverseGstr += 'chr1\t1000\t2000\trow_3\t4\n'
appUniverseGstr += 'chr1\t2000\t3000\trow_4\t12\n'
appUniverseGstr += 'chr2\t0\t1000\trow_6\t7\n'
appUniverseGstr += 'chr2\t10000\t11000\trow_8\t2\n'
appUniverseGstr += 'chr2\t12000\t13000\trow_9\t4\n'

# FUNCTIONS ====================================================================

def test_build_universe_noGroups():
    dom = csd.build([bed1, bed2], 1, csbed)
    with open(dom.fn, "r+") as IH: content = "".join(IH.readlines())
    assert csbed == content

def test_build_universe_groups():
    dom = csd.build([bed1, bed2], 1, csbed, gbed)
    with open(dom.fn, "r+") as IH: content = "".join(IH.readlines())
    assert gstr == content

def test_build_union_noGroups():
    a = bed1.cut(range(3))
    b = bed2.cut(range(3))
    regions = []
    with open(a.fn, "r+") as IH: regions.extend(IH.readlines())
    with open(b.fn, "r+") as IH: regions.extend(IH.readlines())
    regions = list(set(regions))
    regions.sort()

    dom = csd.build([bed1, bed2], 2).sort()
    lines = []
    with open(dom.fn, "r+") as IH: lines.extend(IH.readlines())
    lines.sort()

    assert "".join(regions) == "".join(lines)

def test_build_union_groups():
    dom = csd.build([bed1, bed2], 2, groups = gbed).sort()
    with open(dom.fn, "r+") as IH: content = "".join(IH.readlines())
    assert ngunionstr == content

def test_build_separate():
    dom = csd.build([bed1, bed2], 3)
    assert type(None) == type(dom)

def test_build_intersection_noGroups():
    dom = csd.build([bed1, bed2], 4)
    with open(dom.fn, "r+") as IH: content = "".join(IH.readlines())
    assert ngintersecstr == content

def test_build_intersection_groups():
    dom = csd.build([bed1, bed2], 4, groups = gbed)
    with open(dom.fn, "r+") as IH: content = "".join(IH.readlines())
    assert gintersecstr == content

def test_apply_universe_noGroups():
    dom = csd.build([bed1, bed2], 1, csbed)
    out = csd.apply([bed1, bed2], dom)

    with open(out[0].fn, "r+") as IH: content = IH.read()
    assert content == appUniverseNGstr

def test_apply_universe_groups():
    dom = csd.build([bed1, bed2], 1, csbed, gbed)
    out = csd.apply([bed1, bed2], dom)

    with open(out[0].fn, "r+") as IH: content = IH.read()
    assert content == ""

def test_apply_union_noGroups():
    dom = csd.build([bed1, bed2], 2, csbed)
    out = csd.apply([bed1, bed2], dom)

    with open(out[0].fn, "r+") as IH: content = IH.read()
    assert content == ""

def test_apply_union_groups():
    pass

def test_apply_separate():
    pass

def test_apply_intersection_noGroups():
    pass

def test_apply_intersection_groups():
    pass

# END ==========================================================================

################################################################################
