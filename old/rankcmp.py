#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Date: 20180213
# Project: GPSeq validation by FISH
# Description: calculate difference between rankings obtained with GPSeq,
# 			   either by FISH or seq.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
from joblib import Parallel, delayed
import math
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pp
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import pandas as pd
from tqdm import tqdm
import random
from scipy.stats import norm, ks_2samp
import sys
import time
import warnings

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(description = '''
Calculates distance between rankings obtained with GPSeq, either by FISH or
sequencing. For compatibility with FISH, rank type can be one of: 'chr', 'set',
or 'probe'.

When comparing two sequencing-based rankings, use 'chr' for chromosome-wide
bins, and either set or 'probe' for sub-chromosome bins. When working with non-
-'chr' ranking types, a bed file is required.

Please, note that the script can compare only rankings of the same regions,
which should match those in the provided bed file (if any). If the regions from
the two ranks do not match, the script will work on their intersection. If the
intersection is empty, an error message is displayed.
''', formatter_class = argparse.RawDescriptionHelpFormatter)

# Add mandatory arguments
parser.add_argument('rank1', type = str, nargs = 1,
	help = 'Path to first ranking set.')
parser.add_argument('rank2', type = str, nargs = 1,
	help = 'Path to second ranking set.')

# Add arguments with default value
parser.add_argument('-o', type = str, nargs = 1,
	metavar = 'outdir', default = ['.'],
	help = """Path to output directory, created if missing. Default: '.'""")
parser.add_argument('-t', type = str, nargs = 1,
	metavar = 'type', help = """One of the following: chr, set, probe.
	Default: 'chr'""", choices = ['chr', 'set', 'probe'], default = ['chr'])
parser.add_argument('-d', type = str, nargs = 1,
	metavar = 'distance', choices = ['kt', 'ktw'], default = ['ktw'],
	help = """Wither 'kt' (Kendall tau) or 'ktw' (Kendall tau weighted).
	Default: 'ktw'""")
parser.add_argument('-b', '--bed', type = str, nargs = 1,
	metavar = 'path', default = [None],
	help = """Path to bed file. Needed only for type others than 'chr'.""")
parser.add_argument('-i', type = int, nargs = 1,
	metavar = 'niter', default = [5000],
	help = """Number of iterations to build the random distribution.
	Default: 5000""")
parser.add_argument('-c', '--threads', type = int, nargs = 1,
	metavar = 'nthreads', default = [1],
	help = """Number of threads for parallelization. Default: 1""")
parser.add_argument('-s', type = str, nargs = 1,
	metavar = 'delimiter', default = ['\t'],
	help = """Ranking file field delimiter. Default: TAB""")
parser.add_argument('-p', '--prefix', type = str, nargs = 1,
	metavar = 'text', default = [''],
	help = """Text for output file name prefix. Default: ''""")

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
rank1_path = args.rank1[0]
rank2_path = args.rank2[0]
outdir = args.o[0]
rtype = args.t[0]
dtype = args.d[0]
bed_path = args.bed[0]
if not 'chr' == rtype and type(None) == type(bed_path):
	sys.exit("ERROR: bed path required for non-chromosomal rank types.")
niter = args.i[0]
nthreads = args.threads[0]
sep = args.s[0]
prefix = args.prefix[0]
if not 0 == len(prefix):
	if not "." == prefix[-1]:
		prefix += "."

dlabs = {'kt' : 'Kendall tau distance', 'ktw' : 'Weighted Kendall tau distance'}

# FUNCTIONS ====================================================================

def kendall_tau(r1, r2):
	'''
	Calculate Kendall tau distance between two rankings
	of the same length. Each ranking is a DataFrame of
	two columns: name and value.
	
	Args:
		r1 (pd.DataFrame): first ranking.
		r2 (pd.DataFrame): second ranking.
	
	Returns:
		float: Kendall tau distance.
	'''
	
	# Double-check rankings shapes
	if 4 != r1.shape[1] + r2.shape[1]:
		print("Missing columns in provided rankings.")
		return 1.
	if 0 == r1.shape[0] or 0 == r2.shape[0]:
		print("Provided empty ranking(s).")
		return 1.

	# Reset column labels
	r1.columns = ['name', 'value']
	r2.columns = ['name', 'value']

	# Order rankings based on value
	r1 = r1.sort_values('value')
	r2 = r2.sort_values('value')

	# Check that rankings have the same size
	if r1.shape != r2.shape:
		print("Provided rankings have different length.")
		return 1.
	l = r1.shape[0]

	# Calculate elements indexes
	idx = np.array([range(l),
		[r1['name'].tolist().index(i) for i in r2['name']]]).transpose()

	# Identify all possible couples
	exg = np.array([(x, y) for x in idx[:,0] for y in idx[:,0]])

	# Calculate number of discordant orders
	n = (np.array(idx[exg[:,0],0]) > np.array(idx[exg[:,1],0])).astype('i')
	n += (np.array(idx[exg[:,0],1]) > np.array(idx[exg[:,1],1])).astype('i')
	n = sum(n == 1) / 2.

	# Normalize
	d = n / (l * (l - 1) / 2.)

	# Output
	return d

def kendall_tau_weighted(r1, r2):
	'''
	Calculate Kendall tau weighted distance
	between two rankings of the same length
	Each ranking is a DataFrame of two columns: name and value
	
	Args:
		r1 (pd.DataFrame): first ranking.
		r2 (pd.DataFrame): second ranking.
	
	Returns:
		float: Kendall tau distance.
	'''
	
	# Double-check rankings shapes
	if 4 != r1.shape[1] + r2.shape[1]:
		print("Missing columns in provided rankings.")
		return 1.
	if 0 == r1.shape[0] or 0 == r2.shape[0]:
		print("Provided empty ranking(s).")
		return 1.

	# Reset column labels
	r1.columns = ['name', 'value']
	r2.columns = ['name', 'value']

	# Order rankings based on value
	r1 = r1.sort_values('value')
	r2 = r2.sort_values('value')

	# Check that rankings have the same size
	if r1.shape != r2.shape:
		print("Provided rankings have different length.")
		return 1.
	l = r1.shape[0]

	# Calculate elements indexes
	idx = np.array([range(l),
		[r1['name'].tolist().index(i) for i in r2['name']]]).transpose()

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

def rand_distr(r1, r2, niter, dtype):
	'''
	Generate random distribution of Kendall tau distances
	by permuting one of the provided rankings (permute values).
	Each ranking is a DataFrame of two columns: name and value.

	Args:
		r1 (pd.DataFrame): first ranking.
		r2 (pd.DataFrame): first ranking.
		niter (int): number of iterations.
		dtype (str): either 'kt' or 'ktw'.

	Returns:
		pd.DataFrame: permutated ranking.
	'''

	# Collect distance functions
	dfuns = {'kt' : kendall_tau, 'ktw' : kendall_tau_weighted}

	# Check dtype
	if not dtype in dfuns.keys():
		print("ERRROR: available distances: %s" % (str(dfuns.keys()),))
		return pd.DataFrame([])

	# Make copy of first ranking to avoid shuffling original
	ra = r1.copy()
	rb = r2.copy()

	# Ranking length
	l = ra.shape[0]

	# Empty list to accept distances
	ds = []

	# Iterate shuffling and distance calculation
	if 'kt' == dtype:
		for i in range(niter):
			ra.ix[:, 1] = np.array(ra.iloc[np.random.permutation(l), 1])
			rb.ix[:, 1] = np.array(rb.iloc[np.random.permutation(l), 1])
			ds.append(kendall_tau(ra, rb))
	elif 'ktw' == dtype:
		for i in range(niter):
			ra.ix[:, 1] = np.array(ra.iloc[np.random.permutation(l), 1])
			rb.ix[:, 1] = np.array(rb.iloc[np.random.permutation(l), 1])
			ds.append(kendall_tau_weighted(ra, rb))

	# Output
	return(ds)

def add_region_name(bed, st):
	'''
	Add the 'item' region name column from the bed file.
	
	Args:
		bed (pd.DataFrame): bed file.
		st (pd.DataFrame): sequencing data.
	
	Returns:
		pd.DataFrame: updated sequencing data.
	'''

	# Make copies of input tables
	st = st.copy()

	# Make unique identifier
	bed_uniID = ["%s_%s_%s" % x for x in zip(
		bed['chr'].values, bed['start'].values, bed['end'].values
	)]
	st_uniID = ["%s_%s_%s" % x for x in zip(
		st['chr'].values, st['start'].values, st['end'].values
	)]

	if any([not i in bed_uniID for i in st_uniID]):
		print("ERROR: some ranking table regions cannot be found " +
			"in the provided bed file.")
		return(None)

	# Add region name to st
	idx = [bed_uniID.index(i) for i in st_uniID]
	st = pd.concat([
		st.ix[:, :3].reset_index().drop('index', 1),
		bed.ix[idx, 'name'].reset_index().drop('index', 1),
		st.ix[:, 3:].reset_index().drop('index', 1)
		], axis = 1).rename(columns={'name' : 'item'})

	# Output
	return(st)

def subset_table(st, ft, col = "item"):
	'''
	Subset a sequencing data table based on the FISH data table,
	based on the col region name column.
	
	Args:
		st (pd.DataFrame): sequencing data.
		ft (pd.DataFrame): FISH data.
		col (str): column name for comparison.
	
	Returns:
		pd.DataFrame: updated sequencing data.
	'''

	# Make copies of input tables
	st = st.copy()
	
	# Subset
	idx = [i in ft[col].values for i in st[col].values]
	st = st.ix[idx, :].reset_index().drop('index', 1)

	# Output
	return(st)

def dist_study(r1, r2, niter, dtype):
	'''
	Calculate distance between two rankings of the same length. Each ranking is
	a DataFrame of two columns: name and value. Then, shuffle rankings (iterate)
	and re-calculate distance to obtain a significance value (p-value). Also,
	produce a plot if requested.
	
	Args:
		r1 (pd.DataFrame): first ranking.
		r2 (pd.DataFrame): second ranking.
		niter (int): number of iterations.
		dtype (str): either 'kt' or 'ktw'.
	
	Returns:
		dict: containing information and random distribution data.
	'''
	
	# Collect distance functions
	dfuns = {'kt' : kendall_tau, 'ktw' : kendall_tau_weighted}
	dlabs = {'kt' : 'Kendall tau distance',
		'ktw' : 'Weighted Kendall tau distance'}

	# Prepare output dictionary
	dout = {
		'dist_labs' : dlabs,
		'dist_type' : dtype
	}

	# Check dtype
	if not dtype in dfuns.keys():
		print("ERRROR: available distances: %s" % (str(dfuns.keys()),))
		return

	# Calculate distance
	dout['dist'] = dfuns[dtype](r1, r2)

	# Prepare random distribution with shuffling
	dout['rand_distr'] = rand_distr(r1, r2, niter, dtype)
	dout['mu'], dout['sigma'] = norm.fit(dout['rand_distr'])

	# Calculate p-value
	dout['Z'] = (dout['dist'] - dout['mu']) / dout['sigma']
	dout['pval'] = norm.cdf(dout['dist'], dout['mu'], dout['sigma'])
	if 0.5 < dout['pval']:
		dout['pval'] = 1 - dout['pval']
	dout['pval'] = dout['pval'] / 0.5

	# Calculate goodness of fit
	obs = np.histogram(dout['rand_distr'], bins = 100, density = True)[0]
	exp = np.array(norm.pdf(np.linspace(0, 1, 100),
		loc = dout['mu'], scale = dout['sigma']))
	(dout['ks'], dout['ks_pval']) = ks_2samp(obs, exp)

	# Output
	return(dout)

def pairwise_study_iterFun(r1_type, r2_type, rt1, rt2, dtype, niter):
	'''
	Single pair study for parallelized pairwise study.

	Args:
		r1_type (str): type of first ranking.
		r2_type (str): type of second ranking.
		rt1 (pd.DataFrame): first ranking table.
		rt2 (pd.DataFrame): second ranking table.
		dtype (str): distance type.
		niter (int): number of iterations for random distribution.

	Returns:
		tuple: ranking types and study output.
	'''

	# Extract rankings
	r1 = get_ranking(rt1, r1_type)
	r2 = get_ranking(rt2, r2_type)

	# Subset rankings
	r1 = subset_table(r1, r2, 'name')
	r2 = subset_table(r2, r1, 'name')

	# Perform comparison
	study = dist_study(r1, r2, niter, dtype)

	# Output
	return(r1_type, r2_type, study)

def pairwise_study_plotFun(d, drand, mean, sigma, titles, xlab):
	'''
	Single study plot.

	Args:
		d (float): calculated distance.
		drand (list): random distribution of distances.
		mean (float): random distribution mean.
		sigma (float): random distribution sigma.
		titles (list): plot title/subtitle list.
		xlab (str): X-axis label.
	
	Returns:
		Figure: plot figure.
	'''

	# Prepare empty plot window
	fig, ax = pp.subplots()

	# Plot histogram
	ax.hist(drand, 40, normed = True, color = '#fddbc7')

	# Overlay gaussian
	x = np.linspace(0, 1, 1000)
	ax.plot(x, norm.pdf(x, loc = mean, scale = sigma),
		linestyle = '--', color = '#ef8a62', linewidth = 2)

	# Add significance thresholds
	ax.axvline((norm.ppf(.005) * sigma) + mean,
		color = '#2166ac', linewidth = 1.5)
	ax.axvline((norm.ppf(.025) * sigma) + mean,
		color = '#67a9cf', linewidth = 1.5)
	ax.axvline((norm.ppf(1 - .005) * sigma) + mean,
		color = '#2166ac', linewidth = 1.5)
	ax.axvline((norm.ppf(1 - .025) * sigma) + mean,
		color = '#67a9cf', linewidth = 1.5)

	# Add current distance
	ax.axvline(d, linestyle = ':', color = '#b2182b', linewidth = 2)

	# Layout format
	pp.xlim(0,1)
	pp.xlabel(xlab)
	pp.ylabel('Frequency')
	pp.suptitle(titles[0], fontsize = 11)
	pp.title(titles[1], fontsize = 8)
	pp.subplots_adjust(left = 0.1, right = 0.95, top = 0.85, bottom = 0.1)

	# Output figure
	return(fig)

def pairwise_study(rt1, rt2, niter, nthreads, dtype, outpath, title, rtypes):
	'''
	Study distance between al pairs of rankings.

	Args:
		rt1 (pd.DataFrame): first ranking table.
		rt2 (pd.DataFrame): second ranking table.
		niter (int): number of iterations.
		nthreads (int): number of threads.
		dtype (str): distance type.
		outpath (str): pdf output path.
		title (str): plot title.
		rtypes (list): list of ranking types (S: sequencing, F: FISH.

	Returns:
		tuple: distance and p-value matrices.
	'''

	# Extract columns
	rt1_cols = rt1.columns.tolist()
	rt1_cols = [rt1_cols[i]
		for i in range(rt1_cols.index('item') + 1, len(rt1_cols))]
	rt2_cols = rt2.columns.tolist()
	rt2_cols = [rt2_cols[i]
		for i in range(rt2_cols.index('item') + 1, len(rt2_cols))]

	# Count all pairs for status
	pairs = [(x, y) for x in rt1_cols for y in rt2_cols]

	# Parallelize
	print("Running %d comparisons on %d threads..." % (len(pairs), nthreads))
	out = Parallel(n_jobs = nthreads, verbose = 11)(
		delayed(pairwise_study_iterFun)(
			r1_type, r2_type, rt1, rt2, dtype, niter)
		for (r1_type, r2_type) in pairs)

	# Plot output
	print("Plotting...")
	dmat = np.zeros((len(rt1_cols), len(rt2_cols)))
	pmat = np.zeros((len(rt1_cols), len(rt2_cols)))
	kmat = np.zeros((len(rt1_cols), len(rt2_cols)))

	# Point to plot output pdf file
	outpdf = PdfPages(outpath)
	
	pbar = tqdm(total = len(out))
	for i in range(len(out)):
		(r1_type, r2_type, std) = out[i]

		# Fill matrices
		i = rt1_cols.index(r1_type)
		j = rt2_cols.index(r2_type)
		dmat[i, j] = std['dist']
		pmat[i, j] = std['pval']
		kmat[i, j] = std['ks_pval']

		# Plot title
		titles = [title,
			"R1: $\\bf{%s}$ [%s]; R2: $\\bf{%s}$ [%s]" % (
				r1_type.replace("_", "."), rtypes[0],
				r2_type.replace("_", "."), rtypes[1]) +
			"\nKS: %.6f; KS p.val: %.2E" % (std['ks'], std['ks_pval']) +
			"\nn.perm=%d; Z = %.6f; p.val: %.2E;" % (
				niter, std['Z'], std['pval'])]

		# Prepare plot
		pp.close('all')
		fig = pairwise_study_plotFun(
			std['dist'], std['rand_distr'], std['mu'], std['sigma'],
			titles, std['dist_labs'][std['dist_type']])

		# Save figure
		fig.savefig(outpdf, format = 'pdf')
		pbar.update(1)

	pbar.close()

	# Close pdf pointer
	outpdf.close()
	pp.close('all')

	# Return matrices
	return(dmat, pmat, kmat)

def list_rtypes(t):
	'''
	Extract ranking types from table.

	Args:
		t (pd.DataFrame): ranking table.

	Returns:
		list: list of ranking types.
	'''

	# Check that item table is there
	if not 'item' in t.columns:
		print("ERRROR: provided ranking table lack mandatory 'item' column.")
		return

	# Check that the requested rank type exists
	item_id = t.columns.tolist().index('item')

	# Output
	return(t.columns[(item_id + 1):].tolist())

def get_ranking(t, rtype):
	'''
	Extract ranking from table.

	Args:
		t (pd.DataFrame): ranking table.
		rtype (str): ranking type (column).

	Returns:
		pd.DataFrame: selected ranking, two columns: name and value.
	'''

	# Check that item table is there
	if not 'item' in t.columns:
		print("ERRROR: provided ranking table lack mandatory 'item' column.")
		return

	# Check that the requested rank type exists
	if not rtype in list_rtypes(t):
		print("ERRROR: available rank types: %s" % (
			str(list_rtypes(t)),))
		return

	# Subset and rename columns
	out = t.ix[:, ['item', rtype]].rename(
		columns = {'item' : 'name', rtype : 'value'})

	# Remove NaNs
	out = out.ix[np.logical_not(np.isnan(out['value'].values)), :]

	# Output
	return(out)

def read_ranking_table(path, rtype, sep, bed_path = None):
	# Read input file
	rt = pd.read_csv(path, sep)

	# Check that ranking type is compatible
	rtypes = ['chr', 'set', 'probe']
	itypes = ['chr', 'label', 'probe_label']
	if not rtype in rtypes:
		print("ERROR. Ranking table type must be one of: %s" % (str(rtypes),))
		return(None)

	# Read bed-file if necessary
	if rtype in rtypes[1:]:
		if type(None) == type(bed_path):
			print("ERROR: ranking type %s requires a bed file." % (rtype,))
			return(None)
		else:
			bed = pd.read_csv(bed_path, '\t',
				names = ['chr', 'start', 'end', 'name'])

	if 'item' in rt.columns: # FISH ranking
		rt_type = "FISH"

		# Subset FISH table
		rt = rt.ix[itypes[rtypes.index(rtype)] == rt['type'],
			0:(rt.shape[1] - 1)].reset_index().drop('index', 1)

		# Add coordinates
		if rtype != rtypes[0]:
			idx = [bed['name'].tolist().index(i) for i in rt['item']]
			rt = pd.concat([bed.ix[idx, :3].reset_index().drop('index', 1),
				rt.reset_index().drop('index', 1)], axis = 1)

	else: # Sequencing ranking
		rt_type = "SEQ"

		if 'chr' == rtype: # Rename chromosome column to item
			rt = rt.rename(columns = {'chr' : 'item'})

		# Add item column
		if rtype != rtypes[0]:
			rt = add_region_name(bed, rt)

	# Output
	return((rt, rt_type))

def heatmap_plot(data, ticks, cb_lab, xticks, yticks):
	'''
	Plot heatmap.

	Args:
		data (np.ndarray): matrix for heatmap.
		ticks (list): list of colorbar tick values, used for heatmap vlims.
		cb_lab (str): colorbar label.
		xticks (list): labels for x-axis ticks.
		yticks (list): labels for y-axis ticks.

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
	pp.xticks(range(len(xticks)), xticks, rotation='vertical', fontsize = 6)
	pp.yticks(range(len(yticks)), yticks, fontsize = 6)
	pp.subplots_adjust(left = 0.25, right = 0.95, top = 0.95, bottom = 0.25)

	# Output
	return(fig)

def settings_confirm(s):
	'''
	Show dialog to confirm input settings.

	Args:
		s (str): settings string.
	'''

	# Show settings
	print(s)

	# Loop until proper reply
	confirmed = False
	while not confirmed:
		# Ask for confirmation and read reply
		print("\nRun the analysis?\nYes (y), Abort (a)")
		ans = raw_input()

		# Check that reply is valid
		if not ans in ["y", "a"]:
			print("Invalid answer.")
			pass

		# Continue if confirmed
		if "y" == ans:
			print("")
			return

		# Abort if aborted
		if "a" == ans:
			sys.exit("Aborted.")

# RUN ==========================================================================

# Settings confirmation --------------------------------------------------------

sset = """
     # SETTINGS #

        R1 : %s
        R2 : %s
""" % (rank1_path, rank2_path)
if type(None) != type(bed_path):
	sset += "       Bed : %s\n" % (bed_path,)
sset += """    Outdir : %s

      Type : %s
      Dist : %s

    n.iter : %d
 n.threads : %d

       sep : '%s'
    prefix : '%s'
""" % (outdir, rtype, dlabs[dtype], niter, nthreads, sep, prefix)
settings_confirm(sset)

# Output directory
if not os.path.isdir(outdir):
	os.mkdir(outdir)
if not os.path.sep == outdir[-1]:
	outdir += os.path.sep

# Save settings as txt
fset = open("%s%s%s.settings.txt" % (outdir, prefix, rtype), 'w+')
fset.write(sset)
fset.close()

# Input ------------------------------------------------------------------------

# Read ranking tables
(rt1, rt1_type) = read_ranking_table(rank1_path, rtype, sep, bed_path)
(rt2, rt2_type) = read_ranking_table(rank2_path, rtype, sep, bed_path)

# Extract ranking types
lrt1 = list_rtypes(rt1)
lrt2 = list_rtypes(rt2)

# Subset tables to each other
rt1 = subset_table(rt1, rt2)
rt2 = subset_table(rt2, rt1)

# Compare ----------------------------------------------------------------------

# Run pairwise study
plot_titles = {
	'chr' : "Chromosome-wide comparison",
	'set' : "Probe set comparison",
	'probe' : "Single probe comparison"
}
(dmat, pmat, kmat) = pairwise_study(rt1, rt2, niter, nthreads, dtype,
	'%s%s%s.%s.norms.pdf' % (outdir, prefix, rtype, dtype),
	plot_titles[rtype], [rt1_type, rt2_type])

# Plot -------------------------------------------------------------------------

outpdf = PdfPages('%s%s%s.%s.heatmaps.pdf' % (outdir, prefix, rtype, dtype))

# Plot distance heatmap
fig = heatmap_plot(dmat, [0, 0.5, 1], dlabs[dtype], lrt2, lrt1)
fig.savefig(outpdf, format = 'pdf')
pp.close('all')

# Plot distance p.val heatmap
fig = heatmap_plot(pmat, [0, 0.01, 0.05], "Distance p-value", lrt2, lrt1)
fig.savefig(outpdf, format = 'pdf')
pp.close('all')

# Plot KS p.val heatmap
fig = heatmap_plot(kmat, [0, 0.01, 0.05], "KS p-value", lrt2, lrt1)
fig.savefig(outpdf, format = 'pdf')
pp.close('all')

outpdf.close()

# Output -----------------------------------------------------------------------

pd.DataFrame(dmat, columns = lrt2, index = lrt1).to_csv(
	"%s%s%s.dist_table.tsv" % (outdir, prefix, rtype), "\t")
pd.DataFrame(pmat, columns = lrt2, index = lrt1).to_csv(
	"%s%s%s.pval_table.tsv" % (outdir, prefix, rtype), "\t")
pd.DataFrame(kmat, columns = lrt2, index = lrt1).to_csv(
	"%s%s%s.kspv_table.tsv" % (outdir, prefix, rtype), "\t")

# END ==========================================================================

################################################################################
