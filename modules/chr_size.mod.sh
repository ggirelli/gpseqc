#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module identifies the chromosomes seen in the current exper-
#   -iment and calculates their size (last position cut).
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

echo -e " Retrieving chromosome sizes ..."

# Retrieve chromosome sizes
chrSize=$(cat ${bedfiles[@]} | grep -v 'track' | datamash -sg1 -t$'\t' max 3)

# Sort by chromosomes
chrSizePath="$outdir/"$prefix"chr_size$suffix.tsv"
echo -e "$chrSize" | gawk -f "$awkdir/add_chr_id.awk" | sort -k1,1n | \
    cut -f2,3 > "$chrSizePath" & pid=$!; wait $pid

# END ==========================================================================

################################################################################
