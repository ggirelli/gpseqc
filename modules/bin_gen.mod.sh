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

echo -e " Generating bins ..."

# Generate bins
if $chrWide; then
    cat "$outdir/"$prefix"chr_size$suffix.tsv" | \
        gawk '{ print $1 "\t" 0 "\t" $2 }' \
        > "$outdir/"$prefix"$descr$suffix.bed" & pid=$!
else
    cat "$outdir/"$prefix"chr_size$suffix.tsv" | \
        gawk -v size=$binSize -v step=$binStep -f "$awkdir/mk_bins.awk" \
        > "$outdir/"$prefix"$descr$suffix.bed" & pid=$!
fi

# END ==========================================================================

################################################################################
