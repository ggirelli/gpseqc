#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module combines multiple condition into a single table.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

# Combine read count by cutsite and condition ----------------------------------
echo -e " Combining information ..."

# Combine
comb=""
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    cond_n_reads=$(cat "${bedfiles[$bfi]}" | grep -v "track" | datamash sum 4)
    tmp=$(cat "${bedfiles[$bfi]}" | gawk -v cnr=$cond_n_reads -v bfi=$bfi \
    	-f "$awkdir/add_cnr_bfi.awk")
    comb="$comb$tmp\n"

    # Remove bin_stats
    if [ $notOriginalBed -a false == $debugging ]; then
    	rm "${bedfiles[$bfi]}"; fi
done

# ------------------------------------------------------------------------------
# Columns of $comb:
# 1   2     3   4        5       6      7         8   9
# chr|start|end|condRead|readSum|readMu|readSigma|nCS|condID
comb=$(echo -e "$comb" | gawk -f "$awkdir/add_chr_id.awk" | \
    sort -k1,1n -k3,3n -k10,10n  | cut -f2- )

# Remove first line if empty ---------------------------------------------------
tmp=$(echo -e "$comb" | head -n1)
if [ -z "$tmp" ]; then comb=$(echo -e "$comb" | sed 1d); fi

# END ==========================================================================

################################################################################
