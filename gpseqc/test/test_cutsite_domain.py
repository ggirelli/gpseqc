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

from gpseqc import cutsite_domain as csd

# PARAMS =======================================================================

bedstr1  = 'chr1\t500\t506\tr1\t3\n'
bedstr1 += 'chr1\t600\t606\tr2\t2\n'
bedstr1 += 'chr1\t1000\t1006\tr3\t4\n'
bedstr1 += 'chr1\t2050\t2056\tr4\t4\n'
bedstr1 += 'chr1\t2100\t2106\tr5\t0\n'
bedstr1 += 'chr2\t0\t6\tr7\t5\n'
bedstr1 += 'chr2\t20\t26\tr8\t2\n'
bedstr1 += 'chr2\t10000\t10006\tr9\t2\n'
bedstr1 += 'chr2\t12000\t12006\tr11\t4\n'
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

unionstr  = 'chr1\t1000\t1006\n'
unionstr += 'chr1\t2000\t2006\n'
unionstr += 'chr1\t2050\t2056\n'
unionstr += 'chr1\t2100\t2106\n'
unionstr += 'chr1\t500\t506\n'
unionstr += 'chr1\t600\t606\n'
unionstr += 'chr2\t0\t6\n'
unionstr += 'chr2\t10\t16\n'
unionstr += 'chr2\t10000\t10006\n'
unionstr += 'chr2\t12000\t12006\n'
unionstr += 'chr2\t20\t26\n'

gunionstr  = 'chr1\t0\t1000\n'
gunionstr += 'chr1\t1000\t2000\n'
gunionstr += 'chr1\t2000\t3000\n'
gunionstr += 'chr2\t0\t1000\n'
gunionstr += 'chr2\t10000\t11000\n'
gunionstr += 'chr2\t12000\t13000\n'

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
appUniverseNGstr += 'chr1\t2000\t2006\trow_4\t0\n'
appUniverseNGstr += 'chr1\t2050\t2056\trow_5\t4\n'
appUniverseNGstr += 'chr1\t2100\t2106\trow_6\t0\n'
appUniverseNGstr += 'chr2\t0\t6\trow_7\t5\n'
appUniverseNGstr += 'chr2\t10\t16\trow_8\t0\n'
appUniverseNGstr += 'chr2\t20\t26\trow_9\t2\n'
appUniverseNGstr += 'chr2\t10000\t10006\trow_10\t2\n'
appUniverseNGstr += 'chr2\t12000\t12006\trow_11\t4\n'
appUniverseNGstr += 'chr2\t12100\t12106\trow_12\t0\n'

appUniverseGstr  = 'chr1\t0\t1000\trow_1\t5\n'
appUniverseGstr += 'chr1\t1000\t2000\trow_3\t4\n'
appUniverseGstr += 'chr1\t2000\t3000\trow_4\t4\n'
appUniverseGstr += 'chr2\t0\t1000\trow_6\t7\n'
appUniverseGstr += 'chr2\t1000\t2000\trow_8\t0\n'
appUniverseGstr += 'chr2\t2000\t3000\trow_9\t0\n'
appUniverseGstr += 'chr2\t3000\t4000\trow_10\t0\n'
appUniverseGstr += 'chr2\t4000\t5000\trow_11\t0\n'
appUniverseGstr += 'chr2\t5000\t6000\trow_12\t0\n'
appUniverseGstr += 'chr2\t6000\t7000\trow_13\t0\n'
appUniverseGstr += 'chr2\t7000\t8000\trow_14\t0\n'
appUniverseGstr += 'chr2\t8000\t9000\trow_15\t0\n'
appUniverseGstr += 'chr2\t9000\t10000\trow_16\t0\n'
appUniverseGstr += 'chr2\t10000\t11000\trow_17\t2\n'
appUniverseGstr += 'chr2\t11000\t12000\trow_18\t0\n'
appUniverseGstr += 'chr2\t12000\t13000\trow_19\t4\n'

appUnionNGstr  = 'chr1\t500\t506\trow_1\t3\n'
appUnionNGstr += 'chr1\t600\t606\trow_2\t2\n'
appUnionNGstr += 'chr1\t1000\t1006\trow_3\t4\n'
appUnionNGstr += 'chr1\t2000\t2006\trow_4\t0\n'
appUnionNGstr += 'chr1\t2050\t2056\trow_5\t4\n'
appUnionNGstr += 'chr1\t2100\t2106\trow_6\t0\n'
appUnionNGstr += 'chr2\t0\t6\trow_7\t5\n'
appUnionNGstr += 'chr2\t10\t16\trow_8\t0\n'
appUnionNGstr += 'chr2\t20\t26\trow_9\t2\n'
appUnionNGstr += 'chr2\t10000\t10006\trow_10\t2\n'
appUnionNGstr += 'chr2\t12000\t12006\trow_11\t4\n'

appUnionGstr  = 'chr1\t0\t1000\trow_1\t5\n'
appUnionGstr += 'chr1\t1000\t2000\trow_3\t4\n'
appUnionGstr += 'chr1\t2000\t3000\trow_4\t4\n'
appUnionGstr += 'chr2\t0\t1000\trow_6\t7\n'
appUnionGstr += 'chr2\t10000\t11000\trow_8\t2\n'
appUnionGstr += 'chr2\t12000\t13000\trow_9\t4\n'

appSeparatestr  ='chr1\t500\t506\tr1\t3\n'
appSeparatestr +='chr1\t600\t606\tr2\t2\n'
appSeparatestr +='chr1\t1000\t1006\tr3\t4\n'
appSeparatestr +='chr1\t2050\t2056\tr4\t4\n'
appSeparatestr +='chr2\t0\t6\tr7\t5\n'
appSeparatestr +='chr2\t20\t26\tr8\t2\n'
appSeparatestr +='chr2\t10000\t10006\tr9\t2\n'
appSeparatestr +='chr2\t12000\t12006\tr11\t4\n'

appIntersectNGstr  = 'chr1\t500\t506\trow_1\t3\n'
appIntersectNGstr += 'chr1\t600\t606\trow_2\t2\n'
appIntersectNGstr += 'chr1\t1000\t1006\trow_3\t4\n'
appIntersectNGstr += 'chr1\t2050\t2056\trow_4\t4\n'
appIntersectNGstr += 'chr1\t2100\t2106\trow_5\t0\n'
appIntersectNGstr += 'chr2\t0\t6\trow_6\t5\n'
appIntersectNGstr += 'chr2\t20\t26\trow_7\t2\n'
appIntersectNGstr += 'chr2\t10000\t10006\trow_8\t2\n'
appIntersectNGstr += 'chr2\t12000\t12006\trow_9\t4\n'

appIntersectGstr  = 'chr1\t0\t1000\trow_1\t5\n'
appIntersectGstr += 'chr1\t1000\t2000\trow_3\t4\n'
appIntersectGstr += 'chr1\t2000\t3000\trow_4\t4\n'
appIntersectGstr += 'chr2\t0\t1000\trow_6\t7\n'
appIntersectGstr += 'chr2\t10000\t11000\trow_8\t2\n'
appIntersectGstr += 'chr2\t12000\t13000\trow_9\t4\n'

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
    dom = csd.build([bed1, bed2], 2).sort()
    lines = []
    with open(dom.fn, "r+") as IH: lines.extend(IH.readlines())
    lines.sort()

    assert unionstr == "".join(lines)

def test_build_union_groups():
    dom = csd.build([bed1, bed2], 2, groups = gbed).sort()
    with open(dom.fn, "r+") as IH: content = "".join(IH.readlines())
    assert gunionstr == content

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
    assert content == appUniverseGstr

def test_apply_union_noGroups():
    dom = csd.build([bed1, bed2], 2, csbed)
    out = csd.apply([bed1, bed2], dom)

    with open(out[0].fn, "r+") as IH: content = IH.read()
    assert content == appUnionNGstr

def test_apply_union_groups():
    dom = csd.build([bed1, bed2], 2, csbed, gbed)
    out = csd.apply([bed1, bed2], dom)

    with open(out[0].fn, "r+") as IH: content = IH.read()
    assert content == appUnionGstr

def test_apply_separate():
    dom = csd.build([bed1, bed2], 3)
    out = csd.apply([bed1, bed2], dom)

    with open(out[0].fn, "r+") as IH: content = IH.read()
    assert content == appSeparatestr

def test_apply_intersection_noGroups():
    dom = csd.build([bed1, bed2], 4, csbed)
    out = csd.apply([bed1, bed2], dom)

    with open(out[0].fn, "r+") as IH: content = IH.read()
    assert content == appIntersectNGstr

def test_apply_intersection_groups():
    dom = csd.build([bed1, bed2], 4, csbed, gbed)
    out = csd.apply([bed1, bed2], dom)

    with open(out[0].fn, "r+") as IH: content = IH.read()
    assert content == appIntersectGstr

# END ==========================================================================

################################################################################
