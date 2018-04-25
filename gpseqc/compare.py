# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: rank comparison methods.
'''

# DEPENDENCIES =================================================================

import os
import pandas as pd

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
            assert_msg = "ipath or df are required to construct a RankTable"
            raise AssertError(assert_msg)

        req_cols = ["chrom", "start", "end"]
        for c in req_cols:
            assert c in self._df.columns, "missing column: '%s'" % c

        assert self._df.shape[1] > 3, "no metrics found"
        assert self._df.shape[0] > 0, "empty rank table"

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
        return(self._df.to_string())

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

        # Identify intersection
        iset = set(self.iter_regions()).intersection(set(b.iter_regions()))
        if 0 == len(iset): return(None)

        # Copy current dataframe and subset it
        df = self._df.copy()
        for i in range(self._df.shape[0]):
            if tuple(self._df.iloc[i, :].tolist()[:3]) not in iset:
                df.drop(i)

        return(RankTable(df = df))

    def intersection(self, b):
        '''Return the intersection of two RankTables as a new RankTable.
        i.e. all region intervals that are in both RankTables with the rankings
        from the first (self). Same as self & b.'''
        return(self & b)


class MetricTable():
    '''Instance of a metric table, with 4 columns: chr, start, end, metric.

    Attributes:
        _df (pd.DataFrame): parsed metric table.
        _metric (str): label of the metric in the table.
    '''

    _df = None
    _metric = None

    def __init__(self, df):
        '''Build metric table instance from pd.DataFrame.

        Args:
            df (pd.DataFrame): a metric table extracted from a RankTable.
        '''
        
        self._df = df
        self._df.columns.values[:3] = ["chrom", "start", "end"]
        self._metric = self._df.columns[3]

        assert_msg = "expected 4 columns, got %d" % self._df.shape[1]
        assert self._df.shape[1] == 4, assert_msg

        assert self._df.shape[0] > 0, "empty metric table"

        assert_msg = "unrecognized metric '%s'" % self._df.columns[-1]
        assert_msg += "\nAvailable: %s" % str(list(CMETRICS.keys()))
        assert self._df.columns[-1] in CMETRICS.keys(), assert_msg

        assert_msg = "found duplicate regions"
        assert len(set(zip(self._df.iloc[:, 0], self._df.iloc[:, 1],
            self._df.iloc[:, 2]))) == self._df.shape[0], assert_msg

    def __getitem__(self, i):
        '''Return the i-th region.'''
        if i >= len(self):
            raise IndexError("MetricTable object index out of range")
        return self._df.iloc[i, :]

    def __iter__(self):
        '''Yield one region at a time.'''
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        '''Return number of regions in MetricTable.'''
        return self._df.shape[0]

    def __str__(self):
        '''String representation.'''
        return(self._df.to_string())

    @property
    def metric(self):
        return self._metric

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

        # Identify intersection
        iset = set(self.iter_regions()).intersection(set(b.iter_regions()))
        if 0 == len(iset): return(None)

        # Copy current dataframe and subset it
        df = self._df.copy()
        for i in range(self._df.shape[0]):
            if tuple(self._df.iloc[i, :].tolist()[:3]) not in iset:
                df.drop(i)

        return(MetricTable(df = df))

    def intersection(self, b):
        '''Return the intersection of two MetricTables as a new MetricTable.
        i.e. all region intervals that are in both MetricTables with the
        rankings from the first (self). Same as self & b.'''
        return(self & b)



# END ==========================================================================

################################################################################
