# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: rank comparison methods.
'''

# DEPENDENCIES =================================================================

from joblib import Parallel, delayed
from matplotlib import pyplot as pp
import numpy as np
import os
import pandas as pd
from scipy.stats import norm, shapiro
from tqdm import tqdm

from ggc.args import check_threads
from gpseqc.centrality import CMETRICS

import time

# FUNCTIONS ====================================================================


class RankTable(object):
    '''Instance of a rank table generated with gpseqc_estimate.

    Attributes:
      _df (pdf.DataFrame): parsed rank table.
      _avail_metrics (list): metrics available in the rank table.
    '''

    _df = None
    _avail_metrics = []

    def __init__(self, path = None, sep = "\t", df = None):
        '''Read and parse rank table file.

        Args:
            path (str): path to rank table file.
            sep (str): field delimiter.
        '''

        if type(None) != type(path):
            assert os.path.isfile(path), "file not found: %s" % path

            # Read file
            self._sep = sep
            self._path = path
            self._df = pd.read_csv(path, sep, header = 0)
        elif type(None) != type(df):
            self._df = df
        else:
            self._df = pd.DataFrame(columns = ["chrom", "start", "end"])

        req_cols = ["chrom", "start", "end"]
        for c in req_cols:
            assert c in self._df.columns, "missing column: '%s'" % c

        assert_msg = "found duplicate regions"
        assert len(set(self._all_regions())) == self._df.shape[0], assert_msg

        self._avail_metrics = self._df.columns[3:]
        for c in self._avail_metrics:
            assert_msg = "unrecognized metric '%s'" % c
            assert_msg += "\nAvailable: %s" % str(list(CMETRICS.keys()))
            assert c in CMETRICS.keys(), assert_msg

    def __getitem__(self, i):
        '''Return the i-th metric table.'''
        if i >= len(self):
            raise IndexError("RankTable object index out of range")
        return MetricTable(self._df.iloc[:, [0, 1, 2, i + 3]])

    def __iter__(self):
        '''Yield one metric table at a time.'''
        return (self[i] for i in range(len(self)))

    def __len__(self):
        '''Return number of metrics in the RankTable.'''
        return len(self._avail_metrics)

    def __str__(self):
        '''String representation.'''
        return self._df.to_string()

    def available_metrics(self):
        '''Yields available metrics.'''
        return (am for am in self._avail_metrics)

    def iter_regions(self):
        '''Yields regions.'''
        return (r for r in self._all_regions())

    def _all_regions(self):
        '''Return the list of region intervals.'''
        return(list(zip(
            self._df.iloc[:, 0],
            self._df.iloc[:, 1],
            self._df.iloc[:, 2]
        )))

    def __and__(self, b):
        """Return the intersection of two RankTables as a new RankTable.
        i.e. all region intervals that are in both RankTables with the rankings
        from the first (self)."""

        assert_msg = "RankTable expected, got '%s'" % type(b)
        assert type(RankTable()) == type(b), assert_msg

        # Identify intersection
        iset = set(self.iter_regions()).intersection(set(b.iter_regions()))
        if 0 == len(iset): return(RankTable())

        # Copy current dataframe and subset it
        df = self._df.copy()
        regs = self._all_regions()
        for i in range(self._df.shape[0]):
            if regs[i] not in iset:
                df = df.drop(i)
        df.index = range(df.shape[0])
        del regs

        return RankTable(df = df)

    def intersection(self, b):
        '''Return the intersection of two RankTables as a new RankTable.
        i.e. all region intervals that are in both RankTables with the rankings
        from the first (self). Same as self & b.'''
        return self & b

    def _subset(self, b):
        '''Returns intersection of A with B and B with A.'''
        a = self.intersection(b)
        b = b.intersection(a)
        return(a, b)

    def _sort_as(self, b):
        '''Returns A with the same region order as B.'''
        a_regions = self._all_regions()
        b_regions = b._all_regions()

        rows = []
        for r in b_regions:
            assert r in a_regions, "unmatched region sets."
            rows.append(self._df.iloc[a_regions.index(r), :])

        a = RankTable(df = pd.concat(rows, 1).transpose())

        return(a)

    def shuffle(self):
        '''Shuffles the regions of a RankTable.

        Returns:
            A new RankTable with shuffled regions.
        '''
        df = self._df.copy()
        a = np.array(self._all_regions())
        np.random.shuffle(a)
        df.iloc[:, :3] = a
        return RankTable(df = df)

    def compare(self, b, dfun = None, shuffle = False, skipSubset = False,
        progress = False, threads = 1):
        '''Calculate the distance between all the MetricTables in the two
        RankTables. Distance is defined as dfun.

        Args:
            b (MetricTable): second rank.
            dfun (fun): distance function with index1, index2, mt1, mt2 input.
            skipSubset (bool): if the two ranks are already subsetted.
            progress (bool): show progress bar.
            threads (int): number of threads for parallelization.

        Returns:
            pd.DataFrame: table with distance between a pair of metrics in each
                          cell.
        '''

        dfun = dKT_iter if type(None) == type(dfun) else dfun
        threads = check_threads(threads)

        # Subset and sort
        a = self
        if not skipSubset: a, b = a._subset(b)
        a = a._sort_as(b)

        # Run comparison -------------------------------------------------------
        
        # Prepare generator of index couples
        def pair_gen(a, b):
            for aidx in range(len(a)):
                for bidx in range(len(b)):
                    yield (aidx, bidx)
        pgen = pair_gen(a, b)

        if 1 == threads: # Single-thread
            if progress: pgen = tqdm(pgen, total = len(a)*len(b))
            dtab = [dfun(aidx, bidx, a, b, shuffle) for (aidx, bidx) in pgen]
        else: # Parallelized
            dtab = Parallel(n_jobs = threads, verbose = 11 if progress else 0)(
                delayed(dfun)(aidx, bidx, a, b, shuffle)
                for (aidx, bidx) in pgen)

        # Reshape distance table
        dtab = np.array(dtab)
        dtab = dtab.reshape((len(a), len(b)))

        # Add labels
        dtab = pd.DataFrame(dtab)
        dtab.columns = ["R2_%s" % m for m in b.available_metrics()]
        dtab.index = ["R1_%s" % m for m in a.available_metrics()]

        return(dtab)

    def build_rand_distr(self, b, dfun = None, niter = 1000, skipSubset = False,
        progress = False, threads = 1):
        '''Builds a random distribution by shuffling the two Ranking tables
        and then comparing them. This process needs to be iterated a large enough
        number of times to produce a proper distribution.

        Args:
            b (MetricTable): second rank.
            dfun (fun): distance function with index1, index2, mt1, mt2 input.
            niter (int): number of iterations to build the random distribution.
            skipSubset (bool): if the two ranks are already subsetted.
            progress (bool): show progress bar.
            threads (int): number of threads for parallelization.
        '''
        
        assert niter >= 1, "at least one iteration is required."

        dfun = dKT_iter if type(None) == type(dfun) else dfun
        threads = check_threads(threads)

        # Apply subsetting if needed
        a = self
        if not skipSubset: a, b = a._subset(b)

        # Build random distribution --------------------------------------------

        # Build iteration generator
        igen = (i for i in range(niter))

        if 1 == threads:
            if progress: igen = tqdm(igen, total = niter)

            # Calculate distance table after niter shuffles
            ds = []
            for i in igen:
                ds.append(self.compare(b, dfun, True, skipSubset))

        else: # Parallelized
            ds = Parallel(n_jobs = threads, verbose = 11 if progress else 0)(
                delayed(compareNshuffle)(a = a, b = b,
                    dfun = dfun, skipSubset = skipSubset)
                for i in igen)

        return(ds)

    def test_comparison(self, b, dfun = None, niter = 1000, skipSubset = False,
        progress = False, threads = 1):
        '''Compare two rank tables, providing p-value for significance against
        random distribution (built with niter iterations), and Shapiro-Wilk
        p-value of Gaussian goodnes-of-fit over the random distribution. The
        goodness-of-fit is required for the proper calculation of the
        significance p-value.

        Args:
            b (MetricTable): second rank.
            dfun (fun): distance function with index1, index2, mt1, mt2 input.
            niter (int): number of iterations to build the random distribution.
            skipSubset (bool): if the two ranks are already subsetted.
            progress (bool): show progress bar.
            threads (int): number of threads for parallelization.

        Returns:
            dict: with the following keys:
                dist: DataFrame with distance between metrics.
                random_distribution: niter distance DataFrames after shuffling..
                Z: Z statistic.
                Zpval: P-value calculated from Z.
                SW: Shapiro-Wilk statistic.
                SWpval: p-value calculated from SW.
        '''
        
        assert niter >= 1, "at least one iteration is required."

        dfun = dKT_iter if type(None) == type(dfun) else dfun
        threads = check_threads(threads)

        # Apply subsetting if needed
        a = self
        if not skipSubset: a, b = a._subset(b)

        # Compare --------------------------------------------------------------

        if progress: print("> Calculating distances...")
        dtab = a.compare(b, dfun, skipSubset = True,
            progress = progress, threads = threads)

        if progress: print("> Building random distribution [n:%d]..." % niter)
        rand_distr = a.build_rand_distr(b, dfun, niter, skipSubset = True,
            progress = progress, threads = threads)

        if progress: print("> Calculating p-value(s)...")

        Zpval_df = dtab.copy()
        Zpval_df[Zpval_df >= 0] = np.nan
        Z_df = Zpval_df.copy()
        SWpval_df = Zpval_df.copy()
        SW_df = Zpval_df.copy()
        
        pgen = ((i, j) for i in dtab.index for j in dtab.columns)
        if progress: pgen = tqdm(pgen, total = dtab.shape[0] * dtab.shape[1])
        for (i, j) in pgen:
            Z, Zpval, SW, SWpval = compare2randDistr(
                dtab.loc[i, j], [d.loc[i, j] for d in rand_distr])
            Z_df.loc[i, j] = Z
            Zpval_df.loc[i, j] = Zpval
            SW_df.loc[i, j] = SW
            SWpval_df.loc[i, j] = SWpval

        return({
            "dist" : dtab, "random_distribution" : rand_distr,
            "Z" : Z_df, "Zpval" : Zpval_df,
            "SW" : SW_df, "SWpval" : SWpval_df
        })

    def calc_KendallTau(self, b, *args, **kwargs):
        '''Calculate the Kendall Tau distance between all the MetricTables in
        the two RankTables. Additional parameters are passed to the self.compare
        function.
        '''
        return self.compare(b, dKT_iter, *args, **kwargs)

    def dKT(self, *args, **kwargs):
        '''Alias for calc_KendallTau.'''
        return self.calc_KendallTau(*args, **kwargs)

    def calc_KendallTau_weighted(self, b, *args, **kwargs):
        '''Calculate the weighted Kendall Tau distance between all the
        MetricTables in the two RankTables. Additional parameters are passed to
        the self.compare function.
        '''
        return self.compare(b, dKTw_iter, *args, **kwargs)

    def dKTw(self, *args, **kwargs):
        '''Alias for calc_KendallTau_weighted.'''
        return self.calc_KendallTau_weighted(*args, **kwargs)


def dKT_iter(aidx, bidx, a, b, shuffle = False):
    '''Single MetricTable comparison with Kendall tau distance for RankTable
    compare function. Needs to be outside the class to be pickled for Parallel.

    Args:
        aidx, bidx (int): index of metric to compare.
        a, b (RankTable).
        shuffle (bool): shuffle metrics before comparing them.
    '''
    a = np.array(a._all_regions())
    b = np.array(b._all_regions())

    if shuffle:
        np.random.shuffle(a)
        np.random.shuffle(b)

    return dKT(a, b)

def dKTw_iter(aidx, bidx, a, b, shuffle = False):
    '''Single MetricTable comparison with Kendall tau weighted distance for
    RankTable compare function. Needs to be outside the class to be pickled for
    Parallel.

    Args:
        aidx, bidx (int): index of metric to compare.
        a, b (RankTable).
        shuffle (bool): shuffle metrics before comparing them.
    '''
    a = a._df.iloc[:, aidx + 3].copy().values
    b = b._df.iloc[:, bidx + 3].copy().values

    if shuffle:
        np.random.shuffle(a)
        np.random.shuffle(b)

    return dKTw(a, b)

def compareNshuffle(a, b, dfun, skipSubset, *args, **kwargs):
    return a.compare(b, dfun, True, skipSubset)

def compare2randDistr(d, rand_distr):
    '''Compares a value "d" with a random distribution (expected to fit a
    Gaussian). Also provides goodness-of-fit for the Gaussian.

    Args:
        d (float): value to compare to random distribution.
        rand_distr (np.ndarray): random distribution.

    Returns:
        tuple: Z, Zpvalue, SW, SWpvalue.
    '''

    mu, sigma = norm.fit(rand_distr)

    Z = (d - mu) / sigma
    Zpval = norm.cdf(d, mu, sigma)
    if .5 < Zpval: Zpval = 1 - Zpval
    Zpval *= 2

    W, Wpval = shapiro(rand_distr)

    return((Z, Zpval, W, Wpval))


class MetricTable(object):
    '''Instance of a metric table, with 4 columns: chr, start, end, metric.
    
    Attributes:
      _df (pd.DataFrame): parsed metric table.
      _metric (str): label of the metric in the table.
    '''

    _df = None
    _metric = None

    def __init__(self, df = None):
        '''Build metric table instance from pd.DataFrame.

        Args:
            df (pd.DataFrame): a metric table extracted from a RankTable.
        '''
        
        if type(None) == type(df):
            self._df = pd.DataFrame(columns = [
                "chrom", "start", "end", list(CMETRICS.keys())[0]])
        else:
            self._df = df

        assert_msg = "expected 4 columns, got %d" % self._df.shape[1]
        assert self._df.shape[1] == 4, assert_msg

        self._df.columns.values[:3] = ["chrom", "start", "end"]
        self._metric = self._df.columns[3]

        assert_msg = "unrecognized metric '%s'" % self._df.columns[-1]
        assert_msg += "\nAvailable: %s" % str(list(CMETRICS.keys()))
        assert self._df.columns[-1] in CMETRICS.keys(), assert_msg

        assert_msg = "found duplicate regions"
        assert len(set(zip(self._df.iloc[:, 0], self._df.iloc[:, 1],
            self._df.iloc[:, 2]))) == self._df.shape[0], assert_msg

        # Remove np.nan and np.inf
        self._df = self._df.drop(np.where(np.isnan(self.mcol.astype('f')))[0])
        self._df = self._df.drop(np.where(np.isinf(self.mcol.astype('f')))[0])

        # Sort (rank)
        self._df = self._df.sort_values(self._metric)
        self._df.index = range(self._df.shape[0])

    def __getitem__(self, i):
        '''Return the i-th region.'''
        if type(0) == type(i):
            if i >= len(self):
                raise IndexError("MetricTable object index out of range")
            return self._df.iloc[i, :]
        elif type(slice(0)) == type(i):
            return(self._df.iloc[i, :])

    def __iter__(self):
        '''Yield one region at a time.'''
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        '''Return number of regions in MetricTable.'''
        return self._df.shape[0]

    def __str__(self):
        '''String representation.'''
        return self._df.to_string()

    @property
    def metric(self):
        return self._metric

    @property
    def mcol(self):
        return self._df.iloc[:, 3]

    def iter_regions(self):
        '''Yields regions.'''
        for r in self._all_regions():
            yield r

    def _all_regions(self):
        '''Return the list of region intervals.'''
        return(list(zip(
            self._df.iloc[:, 0],
            self._df.iloc[:, 1],
            self._df.iloc[:, 2]
        )))

    def __and__(self, b):
        """Return the intersection of two MetricTables as a new MetricTable.
        i.e. all region intervals that are in both MetricTables with the
        rankings from the first (self)."""

        assert_msg = "MetricTable expected, got '%s'" % type(b)
        assert type(MetricTable()) == type(b), assert_msg

        # Identify intersection
        iset = set(self.iter_regions()).intersection(set(b.iter_regions()))
        if 0 == len(iset): return(None)

        # Copy current dataframe and subset it
        df = self._df.copy()
        regs = self._all_regions()
        for i in range(self._df.shape[0]):
            if regs[i] not in iset:
                df = df.drop(i)
        df.index = range(df.shape[0])
        del regs

        return MetricTable(df = df)

    def intersection(self, b):
        '''Return the intersection of two MetricTables as a new MetricTable.
        i.e. all region intervals that are in both MetricTables with the
        rankings from the first (self). Same as self & b.'''
        return self & b

    def _subset(self, b):
        '''Returns intersection of A with B and B with A.'''
        a = self.intersection(b)
        b = b.intersection(a)
        return(a, b)

    def shuffle(self):
        '''Shuffles the regions of a MetricTable.

        Returns:
            A new MetricTable with shuffled regions.
        '''
        df = self._df.copy()
        a = np.array(self._all_regions())
        np.random.shuffle(a)
        df.iloc[:, :3] = a
        return MetricTable(df = df)

    def calc_KendallTau(self, b, skipSubset = False, progress = False):
        '''Calculate Kendall tau distance between two MetricTables.
        The distance is calculated only on the intersection between the tables.

        Args:
            b (MetricTable): second rank.
            skipSubset (bool): if the two ranks are already subsetted.
            progress (bool): show progress bar.
        '''

        # Apply subsetting if needed
        a = self
        if not skipSubset: a, b = a._subset(b)

        return calc_KendallTau(a._all_regions(), b._all_regions(), progress)

    def dKT(self, *args, **kwargs):
        '''Alias for calc_KendallTau.'''
        return self.calc_KendallTau(*args, **kwargs)

    def calc_KendallTau_weighted(self, b, skipSubset = False, progress = False):
        '''Calculate Kendall tau distance between two MetricTables.
        The distance is calculated only on the intersection between the tables.

        Args:
            b (MetricTable): .
            skipSubset (bool): .
        '''

        # Apply subsetting if needed
        a = self
        if not skipSubset: a, b = a._subset(b)

        b_ordered = []
        for r in a.iter_regions():
            b_ordered.append(b.mcol[b._all_regions().index(r)])
        
        return calc_KendallTau_weighted(a.mcol, b_ordered, progress)

    def dKTw(self, *args, **kwargs):
        '''Alias for calc_KendallTau_weighted.'''
        return self.calc_KendallTau_weighted(*args, **kwargs)


def calc_KendallTau(a, b, progress = False):
    ''''''

    rsize = len(a) # Store in a variable for simplicity
    assert rsize == len(b), "unmatched rank length."

    # For each pair of ordered elements in a, if their ordered is inverted
    # in b, then count it. Divide the count by the number of couples.
    b_ordered = [a.index(r) for r in b]

    n = 0
    igen = (i for i in range(rsize))
    if progress: igen = tqdm(igen, total = rsize)
    for i in igen:
        for j in range(i + 1, rsize):
            if b_ordered[i] > b_ordered[j]:
                n += 1

    d = n / ((rsize * (rsize - 1)) / 2)

    return d

def dKT(*args, **kwargs):
    '''Alias for calc_KendallTau.'''
    return calc_KendallTau(*args, **kwargs)

def calc_KendallTau_weighted(a, b, progress = False):
    ''''''

    rsize = len(a) # Store in a variable for simplicity
    assert rsize == len(b), "unmatched rank length."

    a_total = 0
    b_total = 0
    a_diffs = 0
    b_diffs = 0

    igen = (i for i in range(rsize))
    if progress: igen = tqdm(igen, total = rsize)

    for i in igen:
        for j in range(i + 1, rsize):
            if 2 == np.isnan([b[i], b[j]]).sum():
                continue

            a_total += np.absolute(a[i] - a[j])
            b_total += np.absolute(b[i] - b[j])

            if b[i] > b[j]:
                a_diffs += np.absolute(a[i] - a[j])
                b_diffs += np.absolute(b[i] - b[j])

    assert_msg  = "if both ranks have constant weight, please use the "
    assert_msg += "standard Kendall tau distance instead."
    assert 0 != a_total + b_total, assert_msg

    distance  = np.sum(a_diffs / a_total) / 2.
    distance += np.sum(b_diffs / b_total) / 2.
    
    return distance

def dKTw(*args, **kwargs):
    '''Alias for calc_KendallTau_weighted.'''
    return calc_KendallTau_weighted(*args, **kwargs)

def plot_comparison(d, rand_distr, title, xlab):
    '''
    Single study plot.

    Args:
        d (float): calculated distance.
        rand_distr (list): random distribution of distances.
        title (str): plot title.
        xlab (str): X-axis label.
    
    Returns:
        Figure: plot figure.
    '''

    # Prepare empty plot window
    fig, ax = pp.subplots()

    # Fit Gaussian
    mu, sigma = norm.fit(rand_distr)

    # Plot histogram
    ax.hist(rand_distr, 40, density = True, color = '#fddbc7')

    # Overlay gaussian
    x = np.linspace(0, 1, 1000)
    ax.plot(x, norm.pdf(x, loc = mu, scale = sigma),
        linestyle = '--', color = '#ef8a62', linewidth = 2)

    # Add significance thresholds
    ax.axvline((norm.ppf(.005) * sigma) + mu,
        color = '#2166ac', linewidth = 1.5)
    ax.axvline((norm.ppf(.025) * sigma) + mu,
        color = '#67a9cf', linewidth = 1.5)
    ax.axvline((norm.ppf(1 - .005) * sigma) + mu,
        color = '#2166ac', linewidth = 1.5)
    ax.axvline((norm.ppf(1 - .025) * sigma) + mu,
        color = '#67a9cf', linewidth = 1.5)

    # Add current distance
    ax.axvline(d, linestyle = ':', color = '#b2182b', linewidth = 2)

    # Layout format
    pp.xlim(0,1)
    pp.xlabel(xlab)
    pp.ylabel('Frequency')
    pp.suptitle("", fontsize = 11)
    pp.title(title, fontsize = 8)
    pp.subplots_adjust(left = 0.1, right = 0.95, top = 0.85, bottom = 0.1)

    # Output figure
    return(fig)

def plot_heatmap(data, ticks, cb_lab, outpath = None):
    '''Plot heatmap.

    Args:
        data (pd.DataFrame): matrix for heatmap.
        ticks (list): list of colorbar tick values, used for heatmap vlims.
        cb_lab (str): colorbar label.
        outpath (str): path to output pdf.

    Returns:
        Figure: heatmap figure canvas.
    '''

    # Create emtpy figure canvas
    fig, ax = pp.subplots()

    # Plot heatmap
    cax = pp.imshow(data, cmap='hot',
        interpolation='nearest', vmin = ticks[0], vmax = ticks[-1])

    # Add colorbar
    cbar = fig.colorbar(cax, ticks = ticks, extend = "both")
    cbar.set_label(cb_lab)
    cbar.ax.tick_params(labelsize = 6)

    # Adjust plot parameters
    pp.xticks(range(len(data.columns)), data.columns,
        rotation='vertical', fontsize = 6)
    pp.yticks(range(len(data.index)), data.index, fontsize = 6)
    pp.subplots_adjust(left = 0.25, right = 0.95, top = 0.95, bottom = 0.25)

    if type(None) != type(outpath):
        if os.path.isdir(os.path.dirname(outpath)):
            outpath = os.path.splitext(outpath)
            if outpath[-1] != ".pdf": outpath[-1] = ".pdf"
            outpath = "".join(outpath)

            fig.savefig(outpath, format = 'pdf')
            pp.close('all')

    # Output
    return(fig)

# END ==========================================================================

################################################################################
