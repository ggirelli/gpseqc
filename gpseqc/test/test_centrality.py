# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: centrality-related method tests.
'''

# DEPENDENCIES =================================================================

import numpy as np
import pandas as pd

from gpseqc import centrality as c

# PARAMS =======================================================================

df1 = pd.DataFrame([
    ['chr1', 0, 249221236, 75860, 4.3781381658683, 46.5264380581169, 17327, 876501, 1],
    ['chr1', 0, 249221236, 102923, 5.97694541231127, 165.119679583602, 17220, 1045062, 2],
    ['chr1', 0, 249221236, 580975, 29.1551663572038, 414.534732953575, 19927, 6625898, 3],
    ['chr1', 0, 249221236, 502659, 25.1820550072642, 434.349758213242, 19961, 5706581, 4],
    ['chr1', 0, 249221236, 564971, 28.6830989490785, 395.57305703688, 19697, 6380685, 5],
    ['chr1', 0, 249221236, 639962, 32.2903274635451, 317.043533799097, 19819, 7348941, 6]
])
df2 = pd.DataFrame([
    ['chr10', 0, 135503768, 37874, 3.77306236302052, 2.68100540922241, 10038, 876501, 1],
    ['chr10', 0, 135503768, 45549, 4.52818371607516, 3.48548407771986, 10059, 1045062, 2],
    ['chr10', 0, 135503768, 293384, 25.1723723723724, 18.2363780110097, 11655, 6625898, 3],
    ['chr10', 0, 135503768, 246839, 21.0829347454732, 14.8824240340175, 11708, 5706581, 4],
    ['chr10', 0, 135503768, 285805, 24.7022471910112, 20.041079768043, 11570, 6380685, 5],
    ['chr10', 0, 135503768, 332791, 28.681461690942, 22.5135593268717, 11603, 7348941, 6]
])
df3 = pd.DataFrame([
    ['chr11', 0, 35503768, 37874, 3.77306236302052, np.nan, 10038, 876501, 1],
    ['chr11', 0, 35503768, 45549, 4.52818371607516, 3.48548407771986, 10059, 1045062, 2],
    ['chr11', 0, 35503768, 293384, 25.1723723723724, 18.2363780110097, 11655, 6625898, 3],
    ['chr11', 0, 35503768, 246839, 21.0829347454732, 14.8824240340175, 11708, 5706581, 4],
    ['chr11', 0, 35503768, 285805, 24.7022471910112, 20.041079768043, 11570, 6380685, 5],
    ['chr11', 0, 35503768, 332791, 28.681461690942, 22.5135593268717, 11603, 7348941, 6]
])
df1.index = [0 for i in range(df1.shape[0])]
df2.index = [1 for i in range(df2.shape[0])]
df3.index = [2 for i in range(df3.shape[0])]
df1.columns = ['chrom', 'start', 'end', 'sum', 'mean', 'std', 'count', 'cond_nreads', 'cond']
df2.columns = ['chrom', 'start', 'end', 'sum', 'mean', 'std', 'count', 'cond_nreads', 'cond']
df3.columns = ['chrom', 'start', 'end', 'sum', 'mean', 'std', 'count', 'cond_nreads', 'cond']

# FUNCTIONS ====================================================================

def test_calcP():
    assert c.calc_p(df1, 0) == 75860 / (876501 * 17327)
    assert c.calc_p(df1, 1) == 102923 / (1045062 * 17220)
    assert c.calc_p(df1, 2) == 580975 / (6625898 * 19927)

def test_calcPC():
    p1 = 75860 / (876501 * 17327)
    assert c.calc_pc(df1, 0) == p1
    p2 = 102923 / (1045062 * 17220) + p1
    assert c.calc_pc(df1, 1) == p2
    p3 = 580975 / (6625898 * 19927) + p2
    assert c.calc_pc(df1, 2) == p3

def test_calcPR():
    p1 = 75860 / (876501 * 17327)
    assert c.calc_pr(df1, 0) == p1
    p2 = (102923 + 75860) / (1045062 * 17220 + 876501 * 17327)
    assert c.calc_pr(df1, 1) == p2
    p3 = (580975 + 102923 + 75860)
    p3 /= (6625898 * 19927 + 1045062 * 17220 + 876501 * 17327)
    assert c.calc_pr(df1, 2) == p3

def test_calcVar():
    v1 = np.power(46.5264380581169, 2)
    assert c.calc_var(df1, 0) == v1
    v2 = np.power(165.119679583602, 2)
    assert c.calc_var(df1, 1) == v2

def test_calcFF():
    v1 = np.power(46.5264380581169, 2) / 4.3781381658683
    assert c.calc_ff(df1, 0) == v1
    v2 = np.power(165.119679583602, 2) / 5.97694541231127
    assert c.calc_ff(df1, 1) == v2

def test_calcCV():
    v1 = 46.5264380581169 / 4.3781381658683
    assert c.calc_cv(df1, 0) == v1
    v2 = 165.119679583602 / 5.97694541231127
    assert c.calc_cv(df1, 1) == v2

def test_est2p():
    v = c.est_2p(df1, c.calc_p, lambda x, y: x / y)
    assert v == c.calc_p(df1, -1) / c.calc_p(df1, 0)

def test_estF():
    v = sum([c.calc_p(df1, i) / c.calc_p(df1, 0)
        for i in range(1, df1.shape[0])])
    assert c.est_f(df1, c.calc_p, lambda x, y: x / y) == v

def test_estG():
    v = sum([c.calc_p(df1, i) / c.calc_p(df1, i - 1)
        for i in range(1, df1.shape[0])])
    assert c.est_g(df1, c.calc_p, lambda x, y: x / y) == v

def test_binEstimate():
    est = c.bin_estimate(df1, ["prob_2p", "var_f", "roc_g"], False)

    # prob_2p
    p2p = (639962 / (7348941 * 19819)) / (75860 / (876501 * 17327))
    assert p2p == est["prob_2p"].values[0]

    # var_f
    vf  = np.log(np.power(165.119679583602, 2) / np.power(46.5264380581169, 2))
    vf += np.log(np.power(414.534732953575, 2) / np.power(46.5264380581169, 2))
    vf += np.log(np.power(434.349758213242, 2) / np.power(46.5264380581169, 2))
    vf += np.log(np.power(395.57305703688, 2) / np.power(46.5264380581169, 2))
    vf += np.log(np.power(317.043533799097, 2) / np.power(46.5264380581169, 2))
    assert vf == est["var_f"].values[0]

    # roc_g
    rg0 = df1['sum'].values[:1].sum() / (
        df1['count'].values[:1] * df1['cond_nreads'].values[:1]).sum()
    rg1 = df1['sum'].values[:2].sum() / (
        df1['count'].values[:2] * df1['cond_nreads'].values[:2]).sum()
    rg2 = df1['sum'].values[:3].sum() / (
        df1['count'].values[:3] * df1['cond_nreads'].values[:3]).sum()
    rg3 = df1['sum'].values[:4].sum() / (
        df1['count'].values[:4] * df1['cond_nreads'].values[:4]).sum()
    rg4 = df1['sum'].values[:5].sum() / (
        df1['count'].values[:5] * df1['cond_nreads'].values[:5]).sum()
    rg5 = df1['sum'].values[:6].sum() / (
        df1['count'].values[:6] * df1['cond_nreads'].values[:6]).sum()
    rg = rg1 / rg0 + rg2 / rg1 + rg3 / rg2 + rg4 / rg3 + rg5 / rg4
    assert rg == est["roc_g"].values[0]

def test_rank():
    est = c.bin_estimate(pd.concat([df1, df2]),
        ["prob_2p", "var_f", "roc_g"], False)
    rank = c.rank(est, ["prob_2p", "var_f", "roc_g"], False)
    erank = ['chr1:0-249221236', 'chr10:0-135503768']
    assert all(erank == rank['prob_2p'].values)
    assert all(erank[::-1] == rank['var_f'].values)
    assert all(erank[::-1] == rank['roc_g'].values)

    est = c.bin_estimate(pd.concat([df1, df2, df3]),
    ["prob_2p", "var_f", "roc_g"], False)
    rank = c.rank(est, ["prob_2p", "var_f", "roc_g"], False)
    erank = ['chr1:0-249221236', 'chr10:0-135503768', 'chr11:0-35503768']

    assert all(erank == rank['prob_2p'].values)
    assert str([erank[1], erank[0], np.nan])==str(rank['var_f'].values.tolist())
    assert all([erank[1], erank[2], erank[0]] == rank['roc_g'].values)

# END ==========================================================================

################################################################################
