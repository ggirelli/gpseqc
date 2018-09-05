# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: centrality estimation methods.
'''

# DEPENDENCIES =================================================================

from joblib import delayed, Parallel
import numpy as np
import pandas as pd
from tqdm import tqdm

from ggc.args import check_threads

from gpseqc.stats import score_outliers

# CONSTANTS ====================================================================

CMETRICS = {
    'prob_2p' : lambda st: est_2p(st, calc_p, ratio),
    'prob_f'  : lambda st: est_f(st, calc_p, ratio),
    'prob_g'  : lambda st: est_g(st, calc_p, ratio),
    'cor_2p'  : lambda st: est_2p(st, calc_pc, ratio),
    'cor_f'   : lambda st: est_f(st, calc_pc, ratio),
    'cor_g'   : lambda st: est_g(st, calc_pc, ratio),
    'roc_2p'  : lambda st: est_2p(st, calc_pr, ratio),
    'roc_f'   : lambda st: est_f(st, calc_pr, ratio),
    'roc_g'   : lambda st: est_g(st, calc_pr, ratio),
    'var_2p'  : lambda st: est_2p(st, calc_var, logratio),
    'var_f'   : lambda st: est_f(st, calc_var, logratio),
    'ff_2p'   : lambda st: est_2p(st, calc_ff, lambda x, y: x - y),
    'ff_f'    : lambda st: est_f(st, calc_ff, lambda x, y: x - y),
    'cv_2p'   : lambda st: est_2p(st, calc_cv, lambda x, y: x - y),
    'cv_f'    : lambda st: est_f(st, calc_cv, lambda x, y: x - y),
    'svar_2p' : lambda st: est_2p(st, calc_var, ratio),
    'svar_f'  : lambda st: est_f(st, calc_var, ratio),
    'svar_g'  : lambda st: est_g(st, calc_var, ratio),
    'scv_2p'  : lambda st: est_2p(st, calc_cv, ratio),
    'scv_f'   : lambda st: est_f(st, calc_cv, ratio),
    'scv_g'   : lambda st: est_g(st, calc_cv, ratio),
    'sff_2p'  : lambda st: est_2p(st, calc_ff, ratio),
    'sff_f'   : lambda st: est_f(st, calc_ff, ratio),
    'sff_g'   : lambda st: est_g(st, calc_ff, ratio)
}

# FUNCTIONS ====================================================================

def calc_p(st, ci):
    '''Calculate restriction probability.

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        ci (int): condition index.

    Returns:
        float
    '''
    assert ci < st.shape[0], "requested condition (index) not found."
    row = st.iloc[ci, :]
    p = (row['cond_nreads'] * row['count'])
    p = row['sum'] / p if 0 != p else np.nan
    return(p)

def calc_pc(st, ci):
    '''Calculate cumulative restriction probability.

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        ci (int): condition index.

    Returns:
        float
    '''
    assert ci < st.shape[0], "requested condition (index) not found."
    if ci == 0: return(calc_p(st, 0))
    else: return(calc_p(st, ci) + calc_pc(st, ci - 1))

def calc_pr(st, ci):
    '''Calculate probability of cumulative restriction.

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        ci (int): condition index.

    Returns:
        float
    '''
    assert ci < st.shape[0], "requested condition (index) not found."
    t = st.iloc[:(ci + 1), :]
    if ci == 0:
        return(calc_p(st, 0))
    else:
        p = np.sum(t.loc[:, "cond_nreads"].values * t.loc[:, "count"].values)
        p = t.loc[:, "sum"].sum() / p if 0 != p else np.nan
        return(p)

def calc_var(st, ci):
    '''Calculate variance.

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        ci (int): condition index.

    Returns:
        np.float64
    '''
    assert ci < st.shape[0], "requested condition (index) not found."
    return(np.power(st["std"].values[ci], 2))

def calc_ff(st, ci):
    '''Calculate Fano Factor.

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        ci (int): condition index.

    Returns:
        np.float64
    '''
    assert ci < st.shape[0], "requested condition (index) not found."
    ff = st["mean"].values[ci]
    ff = np.power(st["std"].values[ci], 2) / ff if 0 != ff else np.nan
    return(ff)

def calc_cv(st, ci):
    '''Calculate coefficient of variation

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        ci (int): condition index.

    Returns:
        np.float64
    '''
    assert ci < st.shape[0], "requested condition (index) not found."
    cv = st["mean"].values[ci]
    cv = st["std"].values[ci] / cv if 0 != cv else np.nan
    return(cv)

def est_2p(st, f1, f2):
    '''Estimates centrality by combining conditions in a two-points fashion.

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        f1 (fun): function for calculating condition-wise metrics.
        f2 (fun): function for putting conditions together.

    Returns:
        Estimated centrality.
    '''
    a = f1(st, 0)
    b = f1(st, st.shape[0] - 1)
    return(f2(b, a))

def est_f(st, f1, f2):
    '''Estimates centrality by combining conditions in a fixed fashion.

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        f1 (fun): function for calculating condition-wise metrics.
        f2 (fun): function for putting conditions together.

    Returns:
        Estimated centrality.
    '''
    out = 0
    a = f1(st, 0)
    for i in range(1, st.shape[0]):
        b = f1(st, i)
        out += f2(b, a)
    return(out)

def est_g(st, f1, f2):
    '''Estimates centrality by combining conditions in a global fashion.

    Args:
        st (pd.DataFrame): bin-based data frame (over conditions).
        f1 (fun): function for calculating condition-wise metrics.
        f2 (fun): function for putting conditions together.

    Returns:
        Estimated centrality.
    '''
    out = 0
    a = f1(st, 0)
    for i in range(1, st.shape[0]):
        b = f1(st, i)
        out += f2(b, a)
        a = b
    return(out)

def ratio(x, y):
    if 0 == y: return(np.nan)
    return(x / y)

def logratio(x, y):
    if np.isnan(y) or np.isnan(x): return(np.nan)
    if y != 0 and x != 0: return(np.log(x / y))
    else: return(np.nan)

def bin_estimate(df, mlist, progress = True):
    '''Estimate centrality for each bin in a condition combined data frame.

    Args:
        df (pd.DataFrame): multi-condition data frame.
        mlist (list): list of metrics to calculate.
        progress (bool): show progress bar.
    '''

    # Build generator
    indexes = list(set(df.index))
    igen = (i for i in indexes)
    if progress: igen = tqdm(igen, total = len(indexes))

    # Iterate over bins
    odf = [bin_estimate_single(i, df, mlist) for i in igen]

    # Assemble output
    odf = pd.concat(odf, axis = 1).transpose()
    columns = ['chrom', 'start', 'end']
    columns.extend(mlist)
    odf.columns = columns[:odf.shape[1]]
    odf.index = range(odf.shape[0])

    return(odf)

def bin_estimate_parallel(df, mlist, threads, progress = True):
    '''Estimate centrality for each bin in a condition combined data frame.

    Args:
        df (pd.DataFrame): multi-condition data frame.
        mlist (list): list of metrics to calculate.
        progress (bool): show progress bar.
        threads (int): number of threads for parallelization.
    '''

    threads = check_threads(threads)

    # Build generator
    verbose = 10 if progress else 0

    # Iterate over bins
    odf =  Parallel(n_jobs = threads, verbose = verbose)(
        delayed(bin_estimate_single)(i, df, mlist)
        for i in list(set(df.index)))

    # Assemble output
    odf = pd.concat(odf, axis = 1).transpose()
    columns = ['chrom', 'start', 'end']
    columns.extend(mlist)
    odf.columns = columns[:odf.shape[1]]
    odf.index = range(odf.shape[0])

    return(odf)

def bin_estimate_single(i, df, mlist):
    '''Estimate centrality for a bin in a condition combined data frame.

    Args:
        i (int): bin id.
        df (pd.DataFrame): multi-condition data frame.
        mlist (list): list of metrics to calculate.
    '''

    st = df.loc[i, :]
    st.index = range(st.shape[0])

    # Prepare output
    orow = st.ix[0, ['chrom', 'start', 'end']]

    # Calculate requested metrics
    for m in mlist:
        assert_msg = "unrecognized metric '%s'." % m
        assert_msg += "\nAvailable: %s" % str(list(CMETRICS.keys()))
        assert m in CMETRICS.keys(), assert_msg
        orow[m] = CMETRICS[m](st)

    return(orow)

def rank(cdf, mlist, progress = True, chrWide = False):
    '''Rank regions based on centrality estimates.

    Args:
        cdf (pd.DataFrame): estimates dataframe.
        mlist (list): list of metrics to rank.
        progress (bool): show progress bar.
        chrWide (bool): use chromosomes as labels (no start/end position).

    Returns:
        pd.DataFrame:.
    '''

    for m in mlist:
        assert m in cdf.columns, "missing column '%s' in  estimate table." % m

    # Prepare region labels (no start/end if chromosome wide)
    if chrWide: labels = [cdf["chrom"].values[i] for i in cdf.index]
    else: labels = ["%s:%s-%s" % tuple(cdf.ix[i, :].tolist()[:3])
            for i in cdf.index]
    
    # Prepare region generator, with progress bar if requested
    gen = (m for m in mlist if m in cdf.columns)
    if progress: gen = tqdm(gen, total = len(mlist))

    # Rank and label regions
    odf = []
    for m in gen:
        l = []
        for i in np.argsort(cdf[m].values.astype("f")):
            x = float(cdf[m].values[i])
            if not np.isnan(x) and not np.isinf(x):
                l.append(labels[i])
            else:
                l.append(np.nan)
        odf.append(l)
    
    # Assemble ranks into a single DataFrame
    odf = pd.DataFrame(odf).transpose()
    odf.columns = [m for m in mlist if m in cdf.columns]

    return(odf)

def rescale(data, method = "IQR", prob = None, lim = 1.5):
    '''Rescale estimated scores with selected method and prob/lim.

    Args:
        data (pd.DataFrame).
        method (str): outlier score type. "Z" calculates normal scores, "t"
            calculates t-Student scores, "chi2" gives chi-squared scores. "IQR"
            considers only values lower than or equal to the first quartile or
            greater than or equal to the third quartile, and sets the  score to
            the ratio between the distance to the closest quartile and the IQR.
            The score for values between the first and third quartiles is set to
            0. "MAD" gives the difference between each value and the median,
            divided by the median absolute deviation.
        prob (float): if set, the corresponding p-values instead of scores are
            given. If set to 1, p-values are returned. otherwise, a logical
            vector is formed, indicating which values are exceeding the
            specified threshold. "IQR" mode does not support probabilities, use
            the lim argument instead.
        lim (float): this value can be set for "IQR" scores, to form a logical
            vector for scores that exceed the limit (after module).

    Returns:
        pd.DataFrame: rescaled estimated scores.
    '''
    estimateNamelist = data.columns[3:len(data.columns)]

    for estimateName in estimateNamelist:
        outliers = score_outliers(data.loc[:, estimateName], method, prob, lim)
        non_outliers = np.logical_not(outliers)
        points = data.loc[non_outliers, estimateName]
        data.loc[:, estimateName] = data.loc[:, estimateName] - np.min(points)
        points = data.loc[non_outliers, estimateName]
        data.loc[:, estimateName] = data.loc[:, estimateName] / np.max(points)
        data.loc[:, estimateName] = 2**data.loc[:, estimateName]

    return(data)

# END ==========================================================================

################################################################################
