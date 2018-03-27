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

source $moddir/functions.mod.sh

# Estimation -------------------------------------------------------------------
echo -e " Estimating centrality ..."

# Prepare paste string
spaste=""; for i in $(seq 1 ${#bedfiles[@]}); do spaste="$spaste -"; done

# Prepare output table
metrics=$(echo -e "$comb" | cut -f1-3 | uniq)

# Probability metrics ----------------------------------------------------------

# Probability metric
echo -e " > Probability ..."
prob_mat=$(echo -e "$comb" | cut -f4,5,8 | \
    gawk -v type="p" -f "$awkdir/pre_process.awk" | paste $spaste)
containsElement "prob_2p" "a ${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    probability_two_points=$(echo -e "$prob_mat" | \
        gawk -v calc="ratio" -v type="2p" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$probability_two_points"))
fi
containsElement "prob_f" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    probability_fixed=$(echo -e "$prob_mat" | \
        gawk -v calc="ratio" -v type="f" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$probability_fixed"))
fi
containsElement "prob_g" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    probability_global=$(echo -e "$prob_mat" | \
        gawk -v calc="ratio" -v type="g" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$probability_global"))
fi

# Cumulative ratio metric
echo -e " > Cumulative ratio ..."
cumrat_mat="$prob_mat"

containsElement "cor_2p" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    cumrat_two_points=$(echo -e "$cumrat_mat" | \
        gawk -v calc="ratio" -v cumrat=1 -v type="2p" \
        -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$cumrat_two_points"))
fi
containsElement "cor_f" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    cumrat_fixed=$(echo -e "$cumrat_mat" | \
        gawk -v calc="ratio" -v cumrat=1 -v type="f" \
        -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$cumrat_fixed"))
fi
containsElement "cor_g" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    cumrat_global=$(echo -e "$cumrat_mat" | \
        gawk -v calc="ratio" -v cumrat=1 -v type="g" \
        -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$cumrat_global"))
fi

# Ratio cumulative metric
echo -e " > Ratio cumulative ..."
ratcum_mat=$(echo -e "$comb" | cut -f4,5,8 | \
    tr '\t' ',' | paste $spaste)

containsElement "roc_2p" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    ratcum_two_points=$(echo -e "$ratcum_mat" | \
        gawk -v calc="ratio" -v ratcum=1 -v type="2p" \
        -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$ratcum_two_points"))
fi
containsElement "roc_f" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    ratcum_fixed=$(echo -e "$ratcum_mat" | \
        gawk -v calc="ratio" -v ratcum=1 -v type="f" \
        -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$ratcum_fixed"))
fi
containsElement "roc_g" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    ratcum_global=$(echo -e "$ratcum_mat" | \
        gawk -v calc="ratio" -v ratcum=1 -v type="g" \
        -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$ratcum_global"))
fi

# Variability metrics ----------------------------------------------------------

# Variance metric
echo -e " > Variance ..."
var_mat=$(echo -e "$comb" | cut -f7 | paste $spaste)

containsElement "var_2p" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    var_two_points=$(echo -e "$var_mat" | \
        gawk -v calc="logratio" -v type="2p" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$var_two_points"))
fi
containsElement "var_f" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    var_fixed=$(echo -e "$var_mat" | \
        gawk -v calc="logratio" -v type="f" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$var_fixed"))
fi

# Fano factor metric
echo -e " > Fano factor ..."
ff_mat=$(echo -e "$comb" | cut -f6,7 | \
    gawk -v type="ff" -f "$awkdir/pre_process.awk" | paste $spaste)

containsElement "ff_2p" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    ff_two_points=$(echo -e "$ff_mat" | \
        gawk -v calc="diff" -v type="2p" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$ff_two_points"))
fi
containsElement "ff_f" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    ff_fixed=$(echo -e "$ff_mat" | \
        gawk -v calc="diff" -v type="f" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$ff_fixed"))
fi

# Coefficient of variation metric
echo -e " > Coefficient of variation ..."
cv_mat=$(echo -e "$comb" | cut -f6,7 | \
    gawk -v type="cv" -f "$awkdir/pre_process.awk" | paste $spaste)

containsElement "cv_2p" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    cv_two_points=$(echo -e "$cv_mat" | \
        gawk -v calc="diff" -v type="2p" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$cv_two_points"))
fi
containsElement "cv_f" "${calc_metrics[@]}"; if [[ 1 == "$?" ]]; then
    cv_fixed=$(echo -e "$cv_mat" | \
        gawk -v calc="diff" -v type="f" -f "$awkdir/estimate_centrality.awk")
    metrics=$(echo -e "$metrics" | paste -d$'\t' - \
        <(echo -e "$cv_fixed"))
fi

# END ==========================================================================

################################################################################
