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

echo -e " Assigning to bins ..."

# Assign bed reads to bins
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Input/output path
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
    outfile=$(echo -e "$infile" | sed 's/csd/intersected/')

    echo -e " > Assigning reads from $infile ..."
    bedtools intersect -a "$outdir/"$prefix"$descr$suffix.bed" \
        -b "${bedfiles[$bfi]}" -wa -wb -loj | cut -f 1-3,8 | \
        sed 's/-1$/0/' | gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
        > "$outdir/$outfile" & pid=$!

    # Remove nzl file
    wait $pid; if [ false == $debugging ]; then rm ${bedfiles[$bfi]}; fi

    # Point to binned bed file instead of original one
    bedfiles[$bfi]="$outdir/$outfile"
done

# Remove bin bed
if [ false == $debugging ]; then
    rm "$outdir/"$prefix"$descr$suffix.bed"
fi

# END ==========================================================================

################################################################################
