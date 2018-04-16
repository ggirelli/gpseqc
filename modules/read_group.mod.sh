#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module groups reads.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

# Generate groups --------------------------------------------------------------
echo -e " Generating groups ..."

generatedGroupsPath="$outdir/"$prefix"groups.$descr$suffix.bed"
cat "$chrSizePath" | gawk -v size=$groupSize -v step=$groupSize \
    -f "$awkdir/mk_bins.awk" > "$generatedGroupsPath" & pid=$!; wait $pid

# Group bed files --------------------------------------------------------------
echo -e " Grouping reads ..."

# Group every bed file
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Set in/output path
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{print $NF}')
    outfile="$outdir/"$prefix"grouped.$descr.$infile$suffix.tsv"

    # Intersect with -loj to keep empty groups
    bedtools intersect -a "$generatedGroupsPath" -b "${bedfiles[$bfi]}" \
        -wa -wb -loj | cut -f 1-3,8 | sed 's/-1$/0/' | \
        gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
        > "$outfile" & pid=$!; wait $pid

    # Point to group bed file instead of original one
    bedfiles[$bfi]="$outfile"
done
notOriginalBed=true

# END ==========================================================================

################################################################################
