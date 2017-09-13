#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module uses datamash to calculate bin statistics.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

echo -e " Calculating bin statistics ..."

# Stats of beds
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Input/output path
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
    outfile=$(echo -e "$infile" | sed 's/intersected/bin_stats/')

    # Calculate statistics
    echo -e " > Calculating for $infile ..."
    cat "${bedfiles[$bfi]}" | datamash -sg1,2,3 sum 5 mean 5 sstdev 5 count 5 \
        | gawk -f "$awkdir/add_chr_id.awk" | sort -k1,1n -k3,3n | cut -f2- \
        > "$outdir/$outfile" & pid=$!

    # Replace (grouped-)cutsite counts
    if [ -n "$csbed" ]; then
        echo 0
    fi

    # Remove binned
    wait $pid; if [ false == $debugging ]; then rm "${bedfiles[$bfi]}"; fi

    # Point to stats bed file instead of original one
    bedfiles[$bfi]="$outdir/$outfile"
done

# END ==========================================================================

################################################################################
