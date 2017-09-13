#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module ranks the regions based on the estimated centrality.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

echo -e " Ranking bins ..."

ranked=""
n_metrics=$(echo -e "$metrics" | awk '{print NF}' | sort -nu | tail -n 1)
for mi in $(seq 4 $n_metrics); do
    if $chrWide; then
        # Rank chromosomes and skip NaNs
        tmp=$(echo -e "$metrics" | cut -f1,$mi | gawk '"nan" != $4' | \
            sort -k2,2n | cut -f1)
    else
        # Rank sub-chromosome regions and skip NaNs
        tmp=$(echo -e "$metrics" | cut -f1,2,3,$mi | gawk '"nan" != $4' | \
            gawk -f "$awkdir/add_chr_id.awk" | \
            gawk '{ print $1"\t"$2"~"$3"~"$4"\t"$5 }' | \
            sort -k1,1n -k3,3n | cut -f2)
    fi
    if [ -z "$ranked" ]; then
        ranked=$tmp
    else
        ranked=$(paste -d$'\t' <(echo -e "$ranked") <(echo -e "$tmp"))
    fi
done

# END ==========================================================================

################################################################################
