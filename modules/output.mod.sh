#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module writes the final output.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

# Output -----------------------------------------------------------------------
echo -e " Writing output ..."

# Add header -------------------------------------------------------------------

# combalization table
header="chr\tstart\tend"
header="$header\tcondRead\treadSum\treadMu\treadSigma\tnCS\tcondID"
comb=$(echo -e "$comb" | \
    gawk -v header="$header" -f "$awkdir/add_header.awk")

# Metrics table
header="chr\tstart\tend"
header="$header\t$(join_by $'\t' ${calc_metrics[@]})"
metrics=$(echo -e "$metrics" | \
    gawk -v header="$header" -f "$awkdir/add_header.awk")

# Ranked table
header=$(echo -e "$metrics" | head -n1 | cut -f4-)
ranked=$(echo -e "$ranked" | \
    gawk -v header="$header" -f "$awkdir/add_header.awk")

# Write ------------------------------------------------------------------------

# Remove bin positions if chromosome wide
if $chrWide; then
    comb=$(echo -e "$comb" | cut -f1,4-)
    metrics=$(echo -e "$metrics" | cut -f1,4-)
fi

# Write
echo -e "$comb" > "$outdir/"$prefix"combined.$descr$suffix.tsv"
echo -e "$metrics" > "$outdir/"$prefix"estimated.$descr$suffix.tsv"
echo -e "$ranked" > "$outdir/"$prefix"ranked.$descr$suffix.tsv"

# END ==========================================================================

################################################################################
