#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module generates and empty bin bed file.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

# !NOTE! Generated bin bed files should contain only 3 columns, otherwise the
# selected columns from the bedtools intersect command will be messed up, as
# that is performed by index using cut -f.

# Generate bins ----------------------------------------------------------------

generatedBinsPath="$outdir/"$prefix"$descr$suffix.bed"
if [ -n "$binBed" ]; then
	cut -f -3 "$binBed" > "$generatedBinsPath"
else
	echo -e " Generating bins ..."
	if $chrWide; then cat "$chrSizePath" | gawk '{ print $1 "\t" 0 "\t" $2 }' \
	        > "$generatedBinsPath" & pid=$!;
	else cat "$chrSizePath" | gawk -v size=$binSize -v step=$binStep \
	    	-f "$awkdir/mk_bins.awk" > "$generatedBinsPath" & pid=$!; fi
	wait $pid
fi

# END ==========================================================================

################################################################################
