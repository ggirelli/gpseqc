# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: rank comparison methods.
'''

# DEPENDENCIES =================================================================

import numpy as np
import os
import pandas as pd
from tqdm import tqdm

from memory_profiler import profile

from gpseqc.centrality import CMETRICS

# FUNCTIONS ====================================================================


class RankTable(object):
    '''Instance of a rank table generated with gpseqc_estimate.

    Attributes:
      _df (pdf.DataFrame): parsed rank table.
      _avail_metrics (list): metrics available in the rank table.
    '''

    _df = None
    _avail_metrics = []

    def __init__(self, ipath = None, sep = "\t", df = None):
        '''Read and parse rank table file.

        Args:
            ipath (str): path to rank table file.
            sep (str): field delimiter.
        '''

        if type(None) != type(ipath):
            assert os.path.isfile(ipath), "file not found: %s" % ipath

            # Read file
            self._sep = sep
            self._path = ipath
            self._df = pd.read_csv(ipath, sep, header = 0)
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
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        '''Return number of metrics in the RankTable.'''
        return len(self._avail_metrics)

    def __str__(self):
        '''String representation.'''
        return self._df.to_string()

    def available_metrics(self):
        '''Yields available metrics.'''
        for am in self._avail_metrics:
            yield am

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
        """Return the intersection of two RankTables as a new RankTable.
        i.e. all region intervals that are in both RankTables with the rankings
        from the first (self)."""

        assert_msg = "RankTable expected, got '%s'" % type(b)
        assert type(RankTable()) == type(b), assert_msg

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

        return RankTable(df = df)

    def intersection(self, b):
        '''Return the intersection of two RankTables as a new RankTable.
        i.e. all region intervals that are in both RankTables with the rankings
        from the first (self). Same as self & b.'''
        return self & b


class MetricTable():
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

    def KendallTau(self, b, skipSubset = False, progress = False):
        '''Calculate Kendall tau distance between two MetricTables.
        The distance is calculated only on the intersection between the tables.

        Args:
            b (MetricTable): second rank.
            skipSubset (bool): if the two ranks are already subsetted.
            progress (bool): show progress bar.
        '''

        # Apply subsetting if needed
        if not skipSubset:
            a = self.intersection(b)
            b = b.intersection(a)

        # Count number of discordant pairs
        n = 0
        bset = set()            # Growing set of regions in B
        regs = a._all_regions() # Regions in A
        igen = range(len(b))
        if progress: igen = tqdm(igen)
        for i in igen:
            # Find the intersection in the sets of elements before the selected
            # one, in both rankings. Elements outside the intersection are
            # part of a discordat pair. Each discordant pair is counted once.
            breg = tuple(b[i].iloc[:3].tolist())
            bset.add(breg)
            aset = set(regs[:(regs.index(breg) + 1)])
            n += len(aset) - len(bset.intersection(aset))

        # Normalize
        d = n / (len(a) * (len(a) - 1) / 2.)

        # Output
        return d

    @profile
    def KendallTau_weighted(self, b, skipSubset = False):
        '''Calculate Kendall tau distance between two MetricTables.
        The distance is calculated only on the intersection between the tables.

        Args:
            b (MetricTable): .
            skipSubset (bool): .
        '''

        # Apply subsetting if needed
        if not skipSubset:
            a = self.intersection(b)
            b = b.intersection(self)

        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Replace with same set-based approach as in KT calculation.
        # Build dictionary to store metrics as values and regions as keys.

        # Calculate elements indexes
        idx = np.array([self._df.index, [self._all_regions().index(x) 
            for x in b.iter_regions()]], dtype = "i").transpose()

        # Identify all possible couples
        exg = np.array([(x, y) for x in idx[:,0] for y in idx[:,0]])

        # Identify discordant orders
        disc1 = np.array(idx[exg[:,0],0]) > np.array(idx[exg[:,1],0])
        disc2 = np.array(idx[exg[:,0],1]) > np.array(idx[exg[:,1],1])
        disc = (disc1.astype('i') + disc2.astype('i')) == 1

        # Calculate weights
        def calc_weight(r, e):
            v1 = r['value'].values[e[:,0]].astype('f')
            v2 = r['value'].values[e[:,1]].astype('f')
            w = abs(v1 - v2) / sum(abs(v1 - v2))
            return(w)
        w1 = calc_weight(r1, exg)
        w2 = calc_weight(r2, exg)

        # Calculate sum of discordant orders weights
        n = sum((w1[disc] + w2[disc]) / 2.)

        # Normalize
        d = n / (sum(w1 + w2) / 2.)

        # Output
        return d

# END ==========================================================================

################################################################################
