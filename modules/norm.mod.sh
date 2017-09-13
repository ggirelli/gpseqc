#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module normalizes read counts based on the last condition.
#   Only groups/locations present in the last condition are kept from the exper-
#   -iment, while the rest is discarded.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================


echo -e " Normalizing over last condition ..."

for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 2")); do
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | \
        gawk '{ print $NF }')
    if [ 0 -ne $groupSize ]; then
        outfile=$(echo -e "$infile" | sed "s/grouped/normlast/")
    else
        outfile=$prefix"normlast.$descr$infile$suffix.tsv"
    fi

    # Intersect and keep only regions present in last condition
    echo -e " > Normalizing $infile ..."
    bedtools intersect -a <(cat "${bedfiles[-1]}" | awk '0 != $5' ) \
        -b "${bedfiles[$bfi]}" -wb | \
        awk 'BEGIN{ OFS = FS = "\t"; } { $10 = $10 / $5; print $0; }' | \
        cut -f6- > "$outdir/$outfile" & pid=$!

    # Remove grouped bed file
    wait $pid; if [ false == $debugging ]; then rm "${bedfiles[$bfi]}"; fi

    # Point to non-zero-loci bed file instead of original one
    bedfiles[$bfi]="$outdir/$outfile"
done

# Remove last condition
if [ false == $debugging ]; then rm "${bedfiles[-1]}"; fi
unset 'bedfiles[-1]'

# END ==========================================================================

################################################################################
