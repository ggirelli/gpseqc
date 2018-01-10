#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module uses bedtools intersect to assign reads to bins.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

# Assign bed reads to bins -----------------------------------------------------
echo -e " Assigning to bins ..."

for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Input/output path
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
    outfile=$(echo -e "$infile" | sed 's/csd/intersected/')

    # Binning reads
    echo -e " > Assigning reads from $infile ..."
    bedtools intersect -a "$outdir/"$prefix"$descr$suffix.bed" \
        -b "${bedfiles[$bfi]}" -wa -wb -loj | cut -f 1-3,8 | \
        sed 's/-1$/0/' | gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
        > "$outdir/$outfile" & pid=$!; wait $pid

    # Remove nzl file
    if [ $notOriginalBed -a false == $debugging ]; then
        rm "${bedfiles[$bfi]}"; fi

    # Point to binned bed file instead of original one
    bedfiles[$bfi]="$outdir/$outfile"
done
notOriginalBed=true

# Remove bin bed
if [ false == $debugging ]; then rm "$generatedBinsPath"; fi

# END ==========================================================================

################################################################################
