#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module prepares the unit domain for the analysis, whether
#    for grouped or single cutsites.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

# Clean
wait $pid; if [ false == $debugging ]; then
    rm "$outdir/"$prefix"chr_size$suffix.tsv"
fi

echo -e " Preparing cutsites ..."

commonDomain=true
case $csMode in
    1) # Universe
        if [ 0 -eq $groupSize ]; then
            echo -e " > Reading cutsite list ..."
            # Read cutsite list
            csbed=$(cat "$csList" | grep -v "track")
        else
            echo -e " > Reading group list ..."
            # Use group list
            csbed=$(cat "$outdir/"$prefix"groups.$descr$suffix.bed" | grep -v "track")
        fi
    ;;
    2) # Union
        echo -e " > Merging cutsite domains ..."

        # Identify cutsites from all conditions
        csbed=$(cat ${bedfiles[@]} | grep -v "track" | cut -f 3 | sort | uniq |\
            gawk -f "$awkdir/add_chr_id.awk" | sort -k1,1n | cut -f2-)
    ;;
    3) # Separated
        # Silence is golden
    ;;
    4) # Intersection
        echo -e " > Intersecting cutsite domains ..."

        # Intersect after removing empty sites/groups
        csbed=$(cat ${bedfiles[@]} | grep -v "track" | gawk '0 != 5' | \
            cut -f 1-3 | sort | uniq -d | gawk -f "$awkdir/add_chr_id.awk" | \
            sort -k1,1n | cut -f2-)
    ;;
    ?)
        msg="!!! ERROR! Unrecognized cutsite domain option: $csMode"
        echo -e "$helps\n$msg"
        exit 1
    ;;
esac

echo -e " > Counting cutsites ..."
if [ -n "$csbed" ]; then
    # Group cutsites
    if [ 0 -ne $groupSize ]; then
        bedtools intersect \
            -a "$outdir/"$prefix"groups.$descr$suffix.bed" -b <(echo -e "$csbed") \
            -wa -c -loj |  cut -f 1-3 | sed 's/-1$/0/' | \
            gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
            > "$outdir/$prefix.cutsite.$descr$suffix.bed" & pid=$!
        wait $pid
    else
        echo -e "$csbed" > "$outdir/$prefix.cutsite.$descr$suffix.bed"
    fi
fi

for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Set input/output paths (csd : CutSite Domain)
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
    if [ false == $normlast ]; then
        if [ 0 -ne $groupSize ]; then
            outfile=$(echo -e "$infile" | sed "s/grouped/csd/")
        else
            outfile=$prefix"csd.$descr.$infile$suffix.tsv"
        fi
    else
        outfile=$(echo -e "$infile" | sed "s/normlast/csd/")
    fi

    if [ -n "$csbed" ]; then
        # Rename current file
        cp "$outdir/$infile" "$outdir/$outfile"
    else
        # Remove zero-loci or empty groups
        cat "${bedfiles[$bfi]}" | gawk '0 != $5' > "$outdir/$outfile" & pid=$!
        wait $pid

        if [ 0 -ne $groupSize ]; then
            bedtools intersect \
                -a "$outdir/"$prefix"groups.$descr$suffix.bed" \
                -b <(cat "$outdir/$outfile" | cut -f1-3) \
                -wa -c -loj | cut -f 1-3 | sed 's/-1$/0/' | \
                gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
                > "$outdir/$prefix.cutsite.$descr$suffix.bed" & pid=$!
            wait $pid;
        fi
    fi

    # Remove grouped bed file
    if [ 0 -ne $groupSize -a false == $debugging ]; then
        rm "$outdir/$infile"
    fi

    # Point to non-zero-loci bed file instead of original one
    bedfiles[$bfi]="$outdir/$outfile"
done

# Clean
if [ false == $debugging -a 0 -ne $groupSize ]; then
    rm "$outdir/"$prefix"groups.$descr$suffix.bed"
fi

# END ==========================================================================

################################################################################
