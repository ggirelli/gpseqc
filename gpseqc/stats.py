# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: stats methods.
'''

# DEPENDENCIES =================================================================

import numpy as np
import pandas as pd
from scipy.stats import chi2, iqr, norm, t
from scipy.stats.mstats import mquantiles

# CONSTANTS ====================================================================

OUTLIER_METHODS = ["Z", "t", "chi2", "IQR", "MAD"]

# FUNCTIONS ====================================================================

def score_outliers(data, method, prob = None, lim = None):
    '''Identify outliers in a 1-dimensional data array.

    Args:
        data (np.ndarray).
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
        np.ndarray: vector of scores, probabilities or logical vector.
    '''

    assert_msg = "unrecognized score type, available: %s" % str(OUTLIER_METHODS)
    assert method in OUTLIER_METHODS, assert_msg

    if type(None) != type(prob):
        assert prob >= 0, "expected probability greater than or equal to 0."
        assert prob <= 1, "expected probability lower than or equal to 1."

    if type(None) != type(lim):
        assert lim >= 0, "expected limit greater than or equal to 0."

    data = data.astype('f') # Force data to float

    n = len(data)
    mean = data.mean()
    std = data.std()
    Z = (data - mean) / std
    
    if method == "Z":
        if type(None) == type(prob): return Z
        pvalues =  norm.cdf(Z)

    if method == "t":
        t_values = Z * np.sqrt(n - 2) / np.sqrt(n - 1 - np.power(Z, 2))
        if type(None) == type(prob): return t_values
        pvalues =  t.cdf(t_values, n - 2)

    if method == "chi2":
        chisquare = np.power(Z, 2)
        if type(None) == type(prob): return chisquare
        pvalues = chi2.cdf(chisquare, 1)

    if method == "IQR":
        data_notNan = data[np.logical_not(np.isnan(data))]
        q1, q3 = mquantiles(data_notNan, [.25, .75])
        IQR = q3 - q1

        IQR_values = np.zeros(n, dtype = np.float64)
        nan_cond = np.isnan(data)
        IQR_values[nan_cond] = lim
        q1_cond = np.logical_and(data <= q1, np.logical_not(nan_cond))
        IQR_values[q1_cond] = (data[q1_cond] - q1) / IQR
        q3_cond = np.logical_and(data >= q3, np.logical_not(nan_cond))
        IQR_values[q3_cond] = (data[q3_cond] - q3) / IQR
        
        if type(None) == type(lim):
            return IQR_values
        else:
            return np.absolute(IQR_values) >= lim

    if method == "MAD":
        median = mquantiles(data, [.5])
        MAD = 1.4826 * np.median(np.absolute(data - median))

        MAD_values = (data - median) / MAD
        if type(None) == type(prob): return MAD_values
        pvalues =  norm.cdf(MAD_values)

    if 1 == prob:
        return pvalues
    else:
        return pvalues >= prob

# END ==========================================================================

################################################################################
