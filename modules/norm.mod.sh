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

# Normalize over last condition ------------------------------------------------
echo -e " Normalizing over last condition ..."

# Get non-empty locations/groups from last condition
lastBed=$(cat "${bedfiles[-1]}" | awk '0 != $5')

# Cycle through the other conditions
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 2")); do
    # Set in/output path
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
    if [ 0 -ne $groupSize ]; then
        outfile="$outdir/"$(echo -e "$infile" | sed "s/grouped/normlast/")
    else
        outfile="$outdir/"$prefix"normlast.$descr$infile$suffix.tsv"
    fi

    # Intersect and keep only regions present in last condition
    echo -e " > Normalizing $infile ..."
    bedtools intersect -a <(echo -e "$lastBed" ) -b "${bedfiles[$bfi]}" -wb | \
        awk 'BEGIN{ OFS = FS = "\t"; } { $10 = $10 / $5; print $0; }' | \
        cut -f6- > "$outfile" & pid=$!; wait $pid

    # Remove grouped bed file
    if [ $notOriginalBed -a false == $debugging ]; then
        rm "${bedfiles[$bfi]}"
    fi

    # Point to non-zero-loci bed file instead of original one
    bedfiles[$bfi]="$outfile"
done
notOriginalBed=true

# Remove last condition
if [ false == $debugging ]; then rm "${bedfiles[-1]}"; fi
unset 'bedfiles[-1]'

# END ==========================================================================

################################################################################
