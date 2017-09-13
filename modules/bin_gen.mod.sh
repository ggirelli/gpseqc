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

# Generate bins ----------------------------------------------------------------
echo -e " Generating bins ..."

generatedBinsPath="$outdir/"$prefix"$descr$suffix.bed"
if $chrWide; then cat "$chrSizePath" | gawk '{ print $1 "\t" 0 "\t" $2 }' \
        > "$generatedBinsPath" & pid=$!;
else cat "$chrSizePath" | gawk -v size=$binSize -v step=$binStep \
    	-f "$awkdir/mk_bins.awk" > "$generatedBinsPath" & pid=$!; fi
wait $pid

# END ==========================================================================

################################################################################
