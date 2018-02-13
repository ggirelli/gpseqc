#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Date: 20170913
# Project: GPSeq - centrality estimation
# Description: estimate genomic region nuclear centrality.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# CONSTANTS ====================================================================

# Code folder
awkdir="`dirname ${BASH_SOURCE}`/awk/"
moddir="`dirname ${BASH_SOURCE}`/modules/"

# PARAMS =======================================================================
source $moddir/input.mod.sh

# RUN ==========================================================================

# 1) Identify chromosome sizes -------------------------------------------------
source $moddir/chr_size.mod.sh

# 2) Generate bins -------------------------------------------------------------
source $moddir/bin_gen.mod.sh

# 3) Group reads ---------------------------------------------------------------
if [ 0 -ne $groupSize ]; then source $moddir/read_group.mod.sh; fi

# 4) Normalize over last conditon ----------------------------------------------
if $normlast; then source $moddir/norm.mod.sh; fi

# 5) Prepare and apply cutsite domain ------------------------------------------
source $moddir/domain.mod.sh

# 6) Intersect with bedtools ---------------------------------------------------
source $moddir/read_bin.mod.sh

# 7) Calculate bin statistics --------------------------------------------------
source $moddir/bin_stat.mod.sh

# 8) Assemble into bin data table ----------------------------------------------
source $moddir/combine.mod.sh

# 9) Estimate centrality -------------------------------------------------------
source $moddir/estimate.mod.sh

# 10) Rank bins ----------------------------------------------------------------
source $moddir/rank.mod.sh

# 11) Output -------------------------------------------------------------------
source $moddir/output.mod.sh

# Final cleanup ----------------------------------------------------------------

if [ false == $debugging ]; then
    rm "$outdir/"$prefix"chr_size$suffix.tsv"
    if [ -e "$generatedCutsitePath" ]; then rm "$generatedCutsitePath"; fi
fi

# END ==========================================================================

################################################################################
