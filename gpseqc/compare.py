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

       	self._avail_metrics = self._df.columns[3:]
        for c in self._avail_metrics:
        	assert_msg = "unrecognized metric '%s'." % c
        	assert_msg += "\nAvailable: %s" % str(list(CMETRICS.keys()))
        	assert c in CMETRICS.keys(), assert_msg

    def __getitem__(self, i):
    	'''Return the i-th metric table.'''
    	raise IndexError("RankTable object index out of range") if i < len(self)
    	return(self._df.iloc[:, [0, 1, 2, i]])

    def __iter__(self):
    	'''Yield one metric table at a time.'''
    	for i in range(len(self)):
    		yield self[i]

    def __len__(self):
    	'''Return number of metrics in the RankTable.'''
    	return(len(self._avail_metrics))

# END ==========================================================================

################################################################################
