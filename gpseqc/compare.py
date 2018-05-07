# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: rank comparison methods.
'''

# DEPENDENCIES =================================================================

from joblib import Parallel, delayed
import numpy as np
import os
import pandas as pd
from tqdm import tqdm

from ggc.args import check_threads
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

    def compare(self, b, dfun = None, skipSubset = False,
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

        # Apply subsetting if needed
        if not skipSubset:
            a = self.intersection(b)
            b = b.intersection(a)
        else:
            a = self

        # Run comparison -------------------------------------------------------
        
        # Prepare generator of index couples
        def pair_gen(a, b):
            for aidx in range(len(a)):
                for bidx in range(len(b)):
                    yield (aidx, bidx)
        pgen = pair_gen(a, b)

        if 1 == threads: # Single-thread
            if progress: pgen = tqdm(pgen, total = len(a)*len(b))
            dtab = [dfun(aidx, bidx, a, b) for (aidx, bidx) in pgen]
        else: # Parallelized
            dtab = Parallel(n_jobs = threads, verbose = 11 if progress else 0)(
                delayed(dfun)(aidx, bidx, a, b)
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
        if not skipSubset:
            a = self.intersection(b)
            b = b.intersection(a)
        else:
            a = self

        # Build random distribution --------------------------------------------

        # Build iteration generator
        igen = (i for i in range(niter))
        if progress: igen = tqdm(igen, total = niter)

        # Calculate distance table after niter shuffles
        ds = []
        for i in igen:
            a = self.shuffle()
            b = b.shuffle()
            ds.append(a.compare(b, dfun, skipSubset, threads = threads))

        return(ds)

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

def dKT_iter(aidx, bidx, a, b):
    return a[aidx].dKT(b[bidx])

def dKTw_iter(aidx, bidx, a, b):
    return a[aidx].dKTw(b[bidx])

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
        if not skipSubset:
            a = self.intersection(b)
            b = b.intersection(a)
        else:
            a = self

        # Count number of discordant pairs -------------------------------------
        n = 0
        bset = set()            # Growing set of regions in B
        regs = a._all_regions() # Regions in A

        # Prepare regions generator
        igen = range(len(b))
        if progress: igen = tqdm(igen)

        # For each region in B, identify the regions with higher rank in A
        # and B and find the intersection between the two sets. The union
        # minus the intersection of the two sets is the number of discordant
        # pairs needed to calculate the Kendall tau distance.
        for i in igen:
            breg = tuple(b[i].iloc[:3].tolist())
            bset.add(breg)
            aset = set(regs[:(regs.index(breg) + 1)])
            n += len(aset) - len(bset.intersection(aset))

        # Normalize
        d = n / (len(a) * (len(a) - 1) / 2.)

        # Output
        return d

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
        if not skipSubset:
            a = self.intersection(b)
            b = b.intersection(self)
        else:
            a = self

        # Count number of discordant pairs -------------------------------------
        n = 0
        w = 0
        d = 0
        bset = set()

        bregs = b._all_regions()
        bdregs = dict([(bregs[i], b[i].iloc[3]) for i in range(len(b))])
        aregs = a._all_regions()
        adregs = dict([(aregs[i], a[i].iloc[3]) for i in range(len(a))])

        # Prepare regions generator
        igen = range(len(b))
        if progress: igen = tqdm(igen)

        def calcLowerBoundWeight(reg, rset, dregs):
            '''Calculate weight when reg has higher weight.

            Args:
              reg (tuple): bed region (chrom, start, end).
              rset (set): set of regions with weight lower than reg's.
              dregs (dict): (region, weight) dictionary.
            '''
            w = np.array([dregs[r] for r in rset])
            nan_cond = w == 0
            if 0 != sum(nan_cond): w[nan_cond] = np.nan
            return np.absolute((w - dregs[reg]) / w)
        
        def calcHigherBoundWeight(reg, rset, dregs):
            '''Calculate weight when reg has lower weight.

            Args:
              reg (tuple): bed region (chrom, start, end).
              rset (set): set of regions with weight higher than reg's.
              dregs (dict): (region, weight) dictionary.
            '''
            w = np.array([dregs[r] for r in rset])
            if dregs[reg] == 0: return np.nan
            return np.absolute((w - dregs[reg]) / dregs[reg])

        def calcTotalWeight(r, aset, bset, aregs, bregs, adregs, bdregs):
            '''Calculate total weight between two sets, for normalization.

            Args:
              r (tuple): (chrom, start, end).
              aset/bset (set): regions set.
              aregs/bregs (list): list of regions.
              adregs/bdregs (dict): (list, estimate) dictionary.

            Returns:
              float
            '''
            w  = np.nansum(calcLowerBoundWeight(r, bset, bdregs))
            w += np.nansum(calcHigherBoundWeight(r, set(bregs) - bset, bdregs))
            w += np.nansum(calcLowerBoundWeight(r, aset, adregs))
            w += np.nansum(calcHigherBoundWeight(r, set(aregs) - aset, adregs))
            w /= 2
            return w

        def calcDiscordantWeight(r, aset, bset, iset, adregs, bdregs):
            '''Calculate weight between discordant pairs in two sets.

            Args:
              r (tuple): (chrom, start, end).
              aset/bset (set): regions set.
              iset (set): intersection between aset and bset.
              adregs/bdregs (dict): (list, estimate) dictionary.

            Returns:
              float
            '''
            wbl = calcLowerBoundWeight(r, bset - iset, bdregs)
            wal = calcLowerBoundWeight(r, aset - iset, adregs)
            return np.sum([np.nansum(x) for x in [wal, wbl]]) / 2

        # For each region in B, identify the regions with higher rank in A
        # and B and find the intersection between the two sets. The union
        # minus the intersection of the two sets is the number of discordant
        # pairs needed to calculate the Kendall tau distance.
        # 
        # In the weighted approach, retain the estimated centrality of each
        # region and use it as weight. As the table is ranked from peripheral to
        # central (increasing estimate value) the current (i-th) metric will
        # always be greater than or equal to any other in the set.
        for i in igen:
            breg = tuple(b[i].iloc[:3].tolist())
            bset.add(breg)

            aidx = aregs.index(breg)
            aset = set(aregs[:(aidx + 1)])
            iset = aset.intersection(bset)

            n += (aidx + 1) - len(iset)
            w += calcTotalWeight(breg, aset, bset, aregs, bregs, adregs, bdregs)
            d += calcDiscordantWeight(breg, aset, bset, iset, adregs, bdregs)

        # Normalize
        d /= w

        # Output
        return d

    def dKTw(self, *args, **kwargs):
        '''Alias for calc_KendallTau_weighted.'''
        return self.calc_KendallTau_weighted(*args, **kwargs)


# END ==========================================================================

################################################################################
