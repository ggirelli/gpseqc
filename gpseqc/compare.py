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
import pyemd
from scipy.stats import percentileofscore, scoreatpercentile
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
    _allow_custom = False

    def __init__(self, path = None, sep = "\t", df = None,
        allow_custom = False):
        '''Read and parse rank table file.

        Args:
            path (str): path to rank table file.
            sep (str): field delimiter.
            df (pd.DataFrame): alternative to path input.
            allow_custom (bool): allow for custom metrics.
        '''

        self._allow_custom = allow_custom

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
        if not allow_custom:
            for c in self._avail_metrics:
                assert_msg = "unrecognized metric '%s'" % c
                assert_msg += "\nAvailable: %s" % str(list(CMETRICS.keys()))
                assert c in CMETRICS.keys(), assert_msg

        # Remove rows with only nan
        keep_ids = []
        for i in range(self._df.shape[0]):
            nan_count = np.isnan(self._df.iloc[i, 3:].values.tolist()).sum()
            if nan_count != self._df.shape[1] - 3:
                keep_ids.append(i)
        self._df = self._df.iloc[keep_ids, :]
        self._df.index = range(self._df.shape[0])

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

        return RankTable(df = df, allow_custom = self._allow_custom)

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

        indexes = []
        for r in b_regions:
            assert r in a_regions, "unmatched region sets."
            indexes.append(a_regions.index(r))

        a = self._df.copy().iloc[indexes, :]
        a = RankTable(df = a, allow_custom = self._allow_custom)

        return(a)

    def shuffle(self):
        '''Shuffles the regions of a RankTable.

        Returns:
            A new RankTable with shuffled regions.
        '''
        df = self._df.copy()
        a = np.array(self._all_regions())
        a = a[np.random.permutation(len(a))]
        df.iloc[:, :3] = a
        return RankTable(df = df, allow_custom = self._allow_custom)

    def compare(self, b, distance = None, shuffle = False, skipSubset = False,
        progress = False, threads = 1):
        '''Calculate the distance between all the MetricTables in the two
        RankTables. Distance is defined as distance.

        Args:
            b (MetricTable): second rank.
            distance (str): distance metric type (see compare.DISTANCE_FUNS).
            skipSubset (bool): if the two ranks are already subsetted.
            progress (bool): show progress bar.
            threads (int): number of threads for parallelization.

        Returns:
            pd.DataFrame: table with distance between a pair of metrics in each
                          cell.
        '''

        distance = "kt" if type(None) == type(distance) else distance
        assert_msg = "unrecognized '%s' distance type. Available: %s" % (
            distance, str(list(DISTANCE_FUNS.keys())))
        assert distance in DISTANCE_FUNS.keys(), assert_msg
        dfun = DISTANCE_FUNS[distance]
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

    def build_rand_distr(self, b, distance = None, niter = 1000,
        skipSubset = False, progress = False, threads = 1):
        '''Builds a random distribution by shuffling the two Ranking tables
        and then comparing them. This process needs to be iterated a large
        enough number of times to produce a proper distribution.

        Args:
            b (MetricTable): second rank.
            distance (str): distance metric type (see compare.DISTANCE_FUNS).
            niter (int): number of iterations to build the random distribution.
            skipSubset (bool): if the two ranks are already subsetted.
            progress (bool): show progress bar.
            threads (int): number of threads for parallelization.
        '''
        
        assert niter >= 1, "at least one iteration is required."

        distance = "kt" if type(None) == type(distance) else distance
        assert_msg = "unrecognized '%s' distance type. Available: %s" % (
            distance, str(list(DISTANCE_FUNS.keys())))
        assert distance in DISTANCE_FUNS.keys(), assert_msg
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
                ds.append(self.compare(b, distance, True, skipSubset))

        else: # Parallelized
            ds = Parallel(n_jobs = threads, verbose = 11 if progress else 0)(
                delayed(compareNshuffle)(a = a, b = b,
                    distance = distance, skipSubset = skipSubset)
                for i in igen)

        return(ds)

    def test_comparison(self, b, distance = None, niter = 1000, skipSubset = False,
        progress = False, threads = 1):
        '''Compare two rank tables, providing p-value for significance against
        random distribution (built with niter iterations), and Shapiro-Wilk
        p-value of Gaussian goodnes-of-fit over the random distribution. The
        goodness-of-fit is required for the proper calculation of the
        significance p-value.

        Args:
            b (MetricTable): second rank.
            distance (str): distance metric type (see compare.DISTANCE_FUNS).
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

        distance = "kt" if type(None) == type(distance) else distance
        assert_msg = "unrecognized '%s' distance type. Available: %s" % (
            distance, str(list(DISTANCE_FUNS.keys())))
        assert distance in DISTANCE_FUNS.keys(), assert_msg
        threads = check_threads(threads)

        # Apply subsetting if needed
        a = self
        if not skipSubset: a, b = a._subset(b)

        # Compare --------------------------------------------------------------

        if progress: print("> Calculating distances...")
        dtab = a.compare(b, distance, skipSubset = True,
            progress = progress, threads = threads)

        if progress: print("> Building random distribution [n:%d]..." % niter)
        rand_distr = a.build_rand_distr(b, distance, niter, skipSubset = True,
            progress = progress, threads = threads)

        if progress: print("> Calculating p-value(s)...")

        pval_df = dtab.copy()
        pval_df[pval_df >= 0] = np.nan
        
        pgen = ((i, j) for i in dtab.index for j in dtab.columns)
        if progress: pgen = tqdm(pgen, total = dtab.shape[0] * dtab.shape[1])
        for (i, j) in pgen:
            data = [d.loc[i, j] for d in rand_distr]

            pval = percentileofscore(data, dtab.loc[i, j], "mean") / 100.
            if pval > 0.5: pval = 1 - pval
            pval *= 2.

            pval_df.loc[i, j] = pval

        return({
            "dist" : dtab,
            "random_distribution" : rand_distr,
            "pval" : pval_df
        })

    def calc_KendallTau(self, b, *args, **kwargs):
        '''Calculate the Kendall Tau distance between all the MetricTables in
        the two RankTables. Additional parameters are passed to the self.compare
        function.
        '''
        return self.compare(b, "kt", *args, **kwargs)

    def calc_KendallTau_weighted(self, b, *args, **kwargs):
        '''Calculate the weighted Kendall Tau distance between all the
        MetricTables in the two RankTables. Additional parameters are passed to
        the self.compare function.
        '''
        return self.compare(b, "ktw", *args, **kwargs)

    def calc_EarthMoversDistance(self, b, *args, **kwargs):
        '''Calculate the Earth Mover's Distance between all the MetricTables in
        the two RankTables. Additional parameters are passed to the self.compare
        function.
        '''
        return self.compare(b, "emd", *args, **kwargs)


def dKT_iter(aidx, bidx, a, b, shuffle = False):
    '''Single MetricTable comparison with Kendall tau distance for RankTable
    compare function. Needs to be outside the class to be pickled for Parallel.

    Args:
        aidx, bidx (int): index of metric to compare.
        a, b (RankTable).
        shuffle (bool): shuffle metrics before comparing them.
    '''
    a_regions = a._all_regions()
    a = [a_regions[i] for i in np.argsort(a._df.iloc[:, aidx + 3].values)]
    b_regions = b._all_regions()
    b = [b_regions[i] for i in np.argsort(b._df.iloc[:, bidx + 3].values)]

    if shuffle:
        a = np.array(a)
        b = np.array(b)
        b = b[np.random.permutation(len(b))]
        a = [tuple(a[i]) for i in range(a.shape[0])]
        b = [tuple(b[i]) for i in range(b.shape[0])]

    return calc_KendallTau(a, b)

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

    a_isorted = np.argsort(a)
    a = a[a_isorted]
    b = b[np.random.permutation(len(b))] if shuffle else b[a_isorted]

    return calc_KendallTau_weighted(a, b)

def dEMD_iter(aidx, bidx, a, b, shuffle = False):
    '''Single MetricTable comparison with Earth Mover's Distance for RankTable
    compare function. Needs to be outside the class to be pickled for Parallel.
    EMD is calculated for items sorted based on the 1st and 2nd rankings
    separately, and then averaged.

    Args:
        aidx, bidx (int): index of metric to compare.
        a, b (RankTable).
        shuffle (bool): shuffle metrics before comparing them.
    '''

    a = a._df.iloc[:, aidx + 3].copy().values.astype(np.float64)
    b = b._df.iloc[:, bidx + 3].copy().values.astype(np.float64)

    nan_condition = np.logical_not(np.logical_or(np.isnan(a), np.isnan(b)))
    a = a[nan_condition]
    b = b[nan_condition]

    a /= a.sum()
    b /= b.sum()

    distance_matrix = mk2DdistanceMatrix(len(a), len(b))

    a_isorted = np.argsort(a)
    b_asorted = b[np.random.permutation(len(b))] if shuffle else b[a_isorted]
    a_asorted = a[a_isorted]

    d  = calc_EarthMoversDistance(a_asorted, b_asorted, distance_matrix)

    b_isorted = np.argsort(b)
    a_bsorted = a[np.random.permutation(len(a))] if shuffle else a[b_isorted]
    b_bsorted = b[b_isorted]

    d += calc_EarthMoversDistance(a_bsorted, b_bsorted, distance_matrix)

    return d / 2.

def mk2DdistanceMatrix(d0, d1):
    ''''''
    distance_matrix = np.zeros((d0, d1))
    for i in range(distance_matrix.shape[0]):
        distance_matrix[i, :] = range(distance_matrix.shape[1])
        distance_matrix[i, :] = np.absolute(distance_matrix[i, :] - i)
    return distance_matrix

def compareNshuffle(a, b, distance, skipSubset, *args, **kwargs):
    return a.compare(b, distance, True, skipSubset)


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
        a = a[np.random.permutation(len(a))]
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

    def calc_EarthMoversDistance(self, b, skipSubset = False):
        '''Calculate Kendall tau distance between two MetricTables.
        The distance is calculated only on the intersection between the tables.

        Args:
            b (MetricTable): .
            skipSubset (bool): .
        '''

        # Apply subsetting if needed
        a = self
        if not skipSubset: a, b = a._subset(b)

        b_asorted = []
        for r in a.iter_regions():
            b_asorted.append(b.mcol[b._all_regions().index(r)])

        a_bsorted = []
        for r in b.iter_regions():
            a_bsorted.append(a.mcol[a._all_regions().index(r)])

        distance_matrix = mk2DdistanceMatrix(len(a), len(b))
        
        d  = calc_EarthMoversDistance(a.mcol, b_asorted, distance_matrix)
        d += calc_EarthMoversDistance(a_bsorted, b.mcol, distance_matrix)

        return d / 2.


def calc_KendallTau(a, b, progress = False):
    '''Calculate Kendall tau distance between two rankings. Each rank is a list
    of sorted elements (ranked). The rank of an element is its index in the
    list. The two list are expected to contain the same elements, possibly in
    different order.

    Args:
        a (list): first rank.
        b (list): second rank.
        progress (bool): show progress bar.
    '''

    rank_size = len(a) # Store in a variable for simplicity
    assert rank_size == len(b), "unmatched rank length."

    # For each element in B, find its index in A.
    b_reindexed = []
    try:
        for r in b:
            b_reindexed.append(a.index(r))
    except ValueError as e:
        print("Element '%s' from rank 2 could not be found in rank 1.")
        raise

    # For each pair of ordered elements in A, if their order is inverted
    # in B, then count a swap. Divide the swap count by the number of couples.
    swap_count = 0
    index_generator = (i for i in range(rank_size))
    if progress: index_generator = tqdm(index_generator, total = rank_size)
    for i in index_generator:
        for j in range(i + 1, rank_size):
            if b_reindexed[i] > b_reindexed[j]:
                swap_count += 1
    d = 2 * swap_count / (rank_size * (rank_size - 1))

    return d

def calc_KendallTau_weighted(a_weights, b_weights, progress = False):
    '''Calculate Kendall tau distance between two rankings. Only the rank
    weights are expected in input. The first weight list should match the first
    ranking order (i.e., be sorted, increasing). The second weight list should
    also match the order of the first ranking (i.e., possibly not sorted).

    Args:
        a_weights (list): weights from 1st ranking, in 1st ranking order.
        b_weights (list): weights from 2nd ranking, in 2nd ranking order.
        progress (bool): show progress bar.
    '''

    rank_size = len(a_weights) # Store in a variable for simplicity
    assert rank_size == len(b_weights), "unmatched rank length."

    a_total_weight = 0
    b_total_weight = 0
    a_swaps_weight = 0
    b_swaps_weight = 0

    igen = (i for i in range(rank_size))
    if progress: igen = tqdm(igen, total = rank_size)

    for i in igen:
        for j in range(i + 1, rank_size):
            if 0 < np.isnan([a_weights[i], a_weights[j]]).sum(): continue
            if 0 < np.isnan([b_weights[i], b_weights[j]]).sum(): continue

            assert a_weights[i] <= a_weights[j], "1st rank is not sorted."

            a_difference = np.absolute(a_weights[i] - a_weights[j])
            b_difference = np.absolute(b_weights[i] - b_weights[j])

            a_total_weight += a_difference
            b_total_weight += b_difference

            if b_weights[i] > b_weights[j]:
                a_swaps_weight += a_difference
                b_swaps_weight += b_difference

    assert_msg  = "if both ranks have constant weight, please use the "
    assert_msg += "standard Kendall tau distance instead."
    assert 0 != a_total_weight, assert_msg
    assert 0 != b_total_weight, assert_msg

    distance  = np.sum(a_swaps_weight / a_total_weight) / 2.
    distance += np.sum(b_swaps_weight / b_total_weight) / 2.
    
    return distance

def calc_EarthMoversDistance(a_weights, b_weights, distance_matrix):
    '''Calculate Earth Mover's Distance between two rankings.

    Args:
        a_weights (np.ndarray): weights from 1st ranking.
        b_weights (np.ndarray): weights from 2nd ranking.
        distance_matrix (np.ndarray): matrix with pair-wise distance between
                                      ranked items.
    '''
    a_weights = np.array(a_weights, dtype = np.float64)
    b_weights = np.array(b_weights, dtype = np.float64)

    d  = pyemd.emd(a_weights, b_weights, distance_matrix,
        extra_mass_penalty = -1.0)
    d /= distance_matrix.max()

    return d


def plot_comparison(d, rand_distr, title, xlab, xlim):
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

    # Plot histogram
    ax.hist(rand_distr, 40, density = True, color = '#fddbc7')

    # Add significance thresholds
    quantiles = scoreatpercentile(rand_distr, [0.5, 2.5, 97.5, 99.5])
    ax.axvline(quantiles[0], color = '#2166ac', linewidth = 1.5)
    ax.axvline(quantiles[1], color = '#67a9cf', linewidth = 1.5)
    ax.axvline(quantiles[2], color = '#67a9cf', linewidth = 1.5)
    ax.axvline(quantiles[3], color = '#2166ac', linewidth = 1.5)

    # Add current distance
    ax.axvline(d, linestyle = ':', color = '#b2182b', linewidth = 2)

    # Layout format
    if type(None) != type(xlim):
        pp.xlim(*xlim)
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

# CONSTANTS ====================================================================

DISTANCE_FUNS = {
    "kt" :  dKT_iter,
    "ktw" : dKTw_iter,
    "emd" : dEMD_iter
}

DISTANCE_LABELS = {
    "kt" :  "Kendall tau distance",
    "ktw" : "Weighted Kendall tau distance",
    "emd" : "Earth Mover's Distance"
}

# END ==========================================================================

################################################################################
