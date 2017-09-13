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

echo -e " Preparing cutsites ..."

# Identify common cutsites/groups lists for domain preparation -----------------
commonDomain=true
case $csMode in
    1) # Universe
        if [ 0 -eq $groupSize ]; then
            echo -e " > Reading cutsite list ..."
            csbed=$(cat "$csList" | grep -v "track")
        else
            echo -e " > Reading group list ..."
            csbed=$(cat "generatedGroupsPath" | grep -v "track")
        fi
    ;;
    2) # Union
        # Identify cutsites from all conditions
        echo -e " > Merging cutsite domains ..."
        csbed=$(cat ${bedfiles[@]} | grep -v "track" | cut -f 3 | sort | uniq \
            | gawk -f "$awkdir/add_chr_id.awk" | sort -k1,1n | cut -f2-)
    ;;
    3) # Separated
        # Silence is golden
    ;;
    4) # Intersection
        # Intersect after removing empty sites/groups
        echo -e " > Intersecting cutsite domains ..."
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

# Group domain if needed and save counts ---------------------------------------
echo -e " > Counting cutsites ..."
generatedCutsitePath="$outdir/$prefix.cutsite.$descr$suffix.bed"
if [ -n "$csbed" ]; then
    # Group cutsites
    if [ 0 -ne $groupSize ]; then
        bedtools intersect \
            -a "generatedGroupsPath" -b <(echo -e "$csbed") \
            -wa -c -loj |  cut -f 1-3 | sed 's/-1$/0/' | \
            gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
            > "$generatedCutsitePath" & pid=$!; wait $pid
    else
        echo -e "$csbed" > "$generatedCutsitePath"
    fi
fi

# Rename or clean bed files ----------------------------------------------------
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Set input/output paths (csd : CutSite Domain)
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
    if [ false == $normlast ]; then
        if [ 0 -ne $groupSize ]; then
            outfile=$(echo -e "$infile" | sed "s/grouped/csd/")
        else outfile=$prefix"csd.$descr.$infile$suffix.tsv"; fi
    else outfile=$(echo -e "$infile" | sed "s/normlast/csd/"); fi

    # Rename current file
    if [ -n "$csbed" ]; then cp "$outdir/$infile" "$outdir/$outfile"
    # Remove zero-loci or empty groups
    else
        cat "${bedfiles[$bfi]}" | gawk '0 != $5' > "$outdir/$outfile" & pid=$!
        wait $pid

        if [ 0 -ne $groupSize ]; then
            bedtools intersect \
                -a "$outdir/"$prefix"groups.$descr$suffix.bed" \
                -b <(cat "$outdir/$outfile" | cut -f1-3) \
                -wa -c -loj | cut -f 1-3 | sed 's/-1$/0/' | \
                gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
                > "$outdir/$prefix.cutsite.$descr$suffix.bed" & pid=$!
            wait $pid
        fi
    fi

    # Remove grouped bed file
    if [ $notOriginalBed -a false == $debugging ]; then
        rm "$outdir/$infile"; fi

    # Point to non-zero-loci bed file instead of original one
    bedfiles[$bfi]="$outdir/$outfile"
    notOriginalBed=true
done

# Clean
if [ false == $debugging -a 0 -ne $groupSize ];then rm "$generatedGroupsPath";fi

# END ==========================================================================

################################################################################
