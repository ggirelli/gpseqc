#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module uses a multi-condition table to estimate regional
#    nuclear centrality.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

echo -e " Estimating centrality ..."

# Prepare paste string
spaste=""; for i in $(seq 1 ${#bedfiles[@]}); do spaste="$spaste -"; done

# Probability metrics ----------------------------------------------------------

# Probability metric
echo -e " > Probability ..."
prob_mat=$(echo -e "$comb" | cut -f4,5,8 | \
    gawk -v type="p" -f "$awkdir/pre_process.awk" | paste $spaste)
probability_two_points=$(echo -e "$prob_mat" | \
    gawk -v calc="ratio" -v type="2p" -f "$awkdir/estimate_centrality.awk")
probability_fixed=$(echo -e "$prob_mat" | \
    gawk -v calc="ratio" -v type="f" -f "$awkdir/estimate_centrality.awk")
probability_global=$(echo -e "$prob_mat" | \
    gawk -v calc="ratio" -v type="g" -f "$awkdir/estimate_centrality.awk")

# Cumulative ratio metric
echo -e " > Cumulative ratio ..."
cumrat_mat="$prob_mat"
cumrat_two_points=$(echo -e "$cumrat_mat" | \
    gawk -v calc="ratio" -v cumrat=1 -v type="2p" \
    -f "$awkdir/estimate_centrality.awk")
cumrat_fixed=$(echo -e "$cumrat_mat" | \
    gawk -v calc="ratio" -v cumrat=1 -v type="f" \
    -f "$awkdir/estimate_centrality.awk")
cumrat_global=$(echo -e "$cumrat_mat" | \
    gawk -v calc="ratio" -v cumrat=1 -v type="g" \
    -f "$awkdir/estimate_centrality.awk")

# Ratio cumulative metric
echo -e " > Ratio cumulative ..."
ratcum_mat=$(echo -e "$comb" | cut -f4,5,8 | \
    tr '\t' ',' | paste $spaste)
ratcum_two_points=$(echo -e "$ratcum_mat" | \
    gawk -v calc="ratio" -v ratcum=1 -v type="2p" \
    -f "$awkdir/estimate_centrality.awk")
ratcum_fixed=$(echo -e "$ratcum_mat" | \
    gawk -v calc="ratio" -v ratcum=1 -v type="f" \
    -f "$awkdir/estimate_centrality.awk")
ratcum_global=$(echo -e "$ratcum_mat" | \
    gawk -v calc="ratio" -v ratcum=1 -v type="g" \
    -f "$awkdir/estimate_centrality.awk")

# Variability metrics ----------------------------------------------------------

# Variance metric
echo -e " > Variance ..."
var_mat=$(echo -e "$comb" | cut -f7 | paste $spaste)
var_two_points=$(echo -e "$var_mat" | \
    gawk -v calc="logratio" -v type="2p" -f "$awkdir/estimate_centrality.awk")
var_fixed=$(echo -e "$var_mat" | \
    gawk -v calc="logratio" -v type="f" -f "$awkdir/estimate_centrality.awk")

# Fano factor metric
echo -e " > Fano factor ..."
ff_mat=$(echo -e "$comb" | cut -f6,7 | \
    gawk -v type="ff" -f "$awkdir/pre_process.awk" | paste $spaste)
ff_two_points=$(echo -e "$ff_mat" | \
    gawk -v calc="diff" -v type="2p" -f "$awkdir/estimate_centrality.awk")
ff_fixed=$(echo -e "$ff_mat" | \
    gawk -v calc="diff" -v type="f" -f "$awkdir/estimate_centrality.awk")

# Coefficient of variation metric
echo -e " > Coefficient of variation ..."
cv_mat=$(echo -e "$comb" | cut -f6,7 | \
    gawk -v type="cv" -f "$awkdir/pre_process.awk" | paste $spaste)
cv_two_points=$(echo -e "$cv_mat" | \
    gawk -v calc="diff" -v type="2p" -f "$awkdir/estimate_centrality.awk")
cv_fixed=$(echo -e "$cv_mat" | \
    gawk -v calc="diff" -v type="f" -f "$awkdir/estimate_centrality.awk")

# Assemble ---------------------------------------------------------------------

# Prepare output table
metrics=$(echo -e "$comb" | cut -f1-3 | uniq | paste -d$'\t' - \
    <(echo -e "$probability_two_points") \
    <(echo -e "$probability_fixed") \
    <(echo -e "$probability_global") \
    <(echo -e "$cumrat_two_points") \
    <(echo -e "$cumrat_fixed") \
    <(echo -e "$cumrat_global") \
    <(echo -e "$ratcum_two_points") \
    <(echo -e "$ratcum_fixed") \
    <(echo -e "$ratcum_global") \
    <(echo -e "$var_two_points") \
    <(echo -e "$var_fixed") \
    <(echo -e "$ff_two_points") \
    <(echo -e "$ff_fixed") \
    <(echo -e "$cv_two_points") \
    <(echo -e "$cv_fixed") \
    )

# END ==========================================================================

################################################################################
