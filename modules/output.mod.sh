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

echo -e " Writing output ..."

# Add header -------------------------------------------------------------------

# combalization table
header="chr\tstart\tend"
header="$header\tcondRead\treadSum\treadMu\treadSigma\tnCS\tcondID"
comb=$(echo -e "$comb" | \
    gawk -v header="$header" -f "$awkdir/add_header.awk")

# Metrics table
header="chr\tstart\tend"
header="$header\tprob_2p\tprob_f\tprob_g"   # Probability
header="$header\tCoR_2p\tCoR_f\tCoR_g"      # Cumulative of Ratio
header="$header\tRoC_2p\tRoC_f\tRoC_g"      # Ratio of Cumulative
header="$header\tvar_2p\tvar_f"             # Variance
header="$header\tff_2p\tff_f"               # Fano factor
header="$header\tcv_2p\tcv_f"               # Coefficient of Variation
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
