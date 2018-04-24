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


class RankTable():
    '''Instance of a rank table generated with gpseqc_estimate.

    Attributes:
      _df (pdf.DataFrame): parsed rank table.
      avail_metrics (list): metrics available in the rank table.
    '''

    _df = None
    _avail_metrics = []

    def __init__(self, ipath, sep = "\t"):
        '''Read and parse rank table file.

        Args:
        	ipath (str): path to rank table file.
        	sep (str): field delimiter.
        '''

        assert os.path.isfile(ipath), "file not found: %s" % ipath

        # Read file
        self._df = pd.read_csv(ipath, sep, header = True)

        req_cols = ["chr", "start", "end"]
        for c in req_cols:
        	assert c in self._df.columns, "missing column: '%s'" % c

        assert self._df.shape[1] > 3, "no metrics found."

        # CHECK THAT NO REGION IS DEFINED MORE THAN ONCE !!!!!!!!!!!!!!!!!!!!!!!

       	self._avail_metrics = self._df.columns[3:]
        for c in self._avail_metrics:
        	assert_msg = "unrecognized metric '%s'." % c
        	assert_msg += "\nAvailable: %s" % str(list(CMETRICS.keys()))
        	assert c in CMETRICS.keys(), assert_msg

    def __getitem__(self, i):
    	'''Return the i-th metric table.'''
    	raise IndexError("RankTable object index out of range") if i < len(self)
    	return(MetricTable(self._df.iloc[:, [0, 1, 2, i]]))

    def __iter__(self):
    	'''Yield one metric table at a time.'''
    	for i in range(len(self)):
    		yield self[i]

    def __len__(self):
    	'''Return number of metrics in the RankTable.'''
    	return(len(self._avail_metrics))


class MetricTable():
	'''Instance of a metric table, with 4 columns: chr, start, end, metric.

	Attributes:
		_df (pd.DataFrame): parsed metric table.
	'''

	_df = None

	def __init__(self, df):
		'''Build metric table instance from pd.DataFrame.

		Args:
			df (pd.DataFrame): a metric table extracted from a RankTable.
		'''
		assert_msg = "expected 4 columns, got %d" % self._df.shape[1]
		assert self._df.shape[1] == 4, assert_msg

    	assert_msg = "unrecognized metric '%s'." % self._df.columns[-1]
    	assert_msg += "\nAvailable: %s" % str(list(CMETRICS.keys()))
		assert self._df.columns[-1] in CMETRICS.keys(), assert_msg

		# CHECK THAT NO REGION IS DEFINED MORE THAN ONCE !!!!!!!!!!!!!!!!!!!!!!!

		self._df.columns[:3] = ["chr", "start", "end"]
		self._df = df


# END ==========================================================================

################################################################################
