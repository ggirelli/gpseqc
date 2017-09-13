#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Date: 20170821
# Project: GPSeq - centrality estimation
# Description: estimate genomic region nuclear centrality.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# PARAMS =======================================================================


# Help string
helps="
usage: ./estimate_centrality.sh [-h][-d][-n]
                                [-c csMode][-l csList]
                                [-s binSize][-p binStep][-g groupSize]
                                [-r prefix][-u suffix]
                                -o outdir [BEDFILE]...

 Description:
  Estimate global centrality. The script performs the following steps:
   (1) Identify & sort chromosomes
   (2) Generate bins
   (3) Group cutsites (intersect)
   (4) Normalize over last condition.
   (5) Prepare domain
   (6) Assign reads to bins (intersect)
   (7) Calculate bin statistics
   (8) Combine condition into a single table
   (9) Estimate centrality
   (10) Rank bins
   (11) Write output
 
 Requirements:
  - bedtools for bin assignment
  - datamash for calculations
  - gawk for text manipulation.

 Notes:
  # Cutsite domain
   The cutsite domain can be specified as follows:
    1 - all genomic cutsites (universe)
    2 - all cutsites restricted in the experiment (union)
    3 - all cutsites restricted in a condition (separate)
    4 - all cutsites restricted in all conditions (intersection)
   Default is 3 (separate). Also, note that if option 1 is selected, an
   additional argument -l is required.
  
   Statistics (mean, variance) metrics take into account only the cutsites
   included in the specified cutsite domain. The same cutsite domain is used for
   all estimates.

   Options 3 and 4 include an empty-cutsites/grups removal step. In this case,
   they are removed before bin assignment, while empty bins are kept. Also,
   normalization is performed after empty bin removal but before bin assignment,
   i.e., either on the grouped or single cutsites.

  # Resolution
   Depending on the sequencing resolution, it might not be feasible to go for
   single-cutsite resolution. Thus, cutsite can be grouped for the statistics
   calculation using the -g option.
 
   In case of sub-chromosome bins, the ranking is done in an ordered
   chromosome-wise manner.

 Mandatory arguments:
  -o outdir     Output folder.
  BEDFILE       At least two (2) GPSeq condition bedfiles, space-separated and
                in increasing order of restriction conditions intensity.
                Expected to be ordered per condition. As BEDFILE is a positional
                argument, it should be provided after any other argument.

 Optional arguments:
  -h            Show this help page.
  -d            Debugging mode: save intermediate results.
  -n            Use last condition for normalization.
  -c csMode     Custite mode (see Notes).
  -l csList     Path to cutsite bed file. Required for -c1 when -g is not used.
  -s binSize    Bin size in bp. Default to chromosome-wide bins.
  -p binStep    Bin step in bp. Default to bin sizeinStep.
  -g groupSize  Group size in bp. Used to group bins for statistics calculation.
                binSize must be divisible by groupSize. Not used by default.
  -r prefix     Output name prefix.
  -u suffix     Output name suffix.
"

# Default values
csMode=3
binSize=0
binStep=0
groupSize=0
chrWide=true
debugging=false
normlast=false

# Parse options
while getopts hdnc:l:s:p:g:o:r:u: opt; do
    case $opt in
        h)
            # Help page
            echo -e "$helps"
            exit 0
        ;;
        d)
            # Debugging mode
            debugging=true
        ;;
        n)
            # Normalize with last condition
            normlast=true
        ;;
        c)
            # Cutsite domain
            if [ $OPTARG -le 0 -o $OPTARG -ge 5 ]; then
                msg="!!! ERROR! Invalid -c option, valid values: 1-4."
                echo -e "$help\n$msg"
                exit 1
            else
                csMode=$OPTARG
            fi
        ;;
        l)
            # Cutsite bed file
            if [ -e $OPTARG -a -n "$OPTARG" ]; then
                csList="$OPTARG"
            else
                msg="!!!ERROR! Invalid -l option, file not found."
                msg="$msg\n    File: $bf"
                echo -e " $helps\n$msg"
                exit 1
            fi
        ;;
        s)
            # Bin size
            if [ $OPTARG -le 0 ]; then
                msg="!!! ERROR! Invalid -s option. Bin size must be > 0."
                echo -e "$help\n$msg"
                exit 1
            else
                binSize=$OPTARG
                chrWide=false
            fi
        ;;
        p)
            # Bin step
            if [ $OPTARG -le 0 ]; then
                msg="!!! ERROR! Invalid -p option. Bin step must be > 0."
                echo -e "$help\n$msg"
                exit 1
            else
                binStep=$OPTARG
            fi
        ;;
        g)
            # Group size
            if [ $OPTARG -le 0 ]; then
                msg="!!! ERROR! Invalid -g option. Group size must be > 0."
                echo -e "$help\n$msg"
                exit 1
            else
                groupSize=$OPTARG
            fi
        ;;
        o)
            # Output directory
            if [ ! -d $OPTARG ]; then
                mkdir -p $OPTARG
            fi
            outdir=$OPTARG
        ;;
        r)
            # Prefix
            prefix=$OPTARG

            # Add trailing dot
            prefix=$(echo -e "$prefix" | sed -r 's/([^\.])$/\1\./' | \
                tr ' ' '_')
        ;;
        u)
            # Suffix
            suffix=$OPTARG

            # Add leading dot
            suffix=$(echo -e "$suffix" | sed -r 's/^([^\.])/\.\1/' | tr ' ' '_')
        ;;
        ?)
            msg="!!! ERROR! Unrecognized option."
            echo -e "$help\n$msg"
            exit 1
        ;;
    esac
done

# Check mandatory options
if [ -z "$outdir" ]; then
    echo -e "$helps\n!!! ERROR! Missing mandatory -o option.\n"
    exit 1
fi
if [ ! -x "$(command -v bedtools)" -o -z "$(command -v bedtools)" ]; then
    echo -e "$helps\n!!! ERROR! Missing bedtools.\n"
    exit 1
fi
if [ ! -x "$(command -v datamash)" -o -z "$(command -v datamash)" ]; then
    echo -e "$helps\n!!! ERROR! Missing datamash.\n"
    exit 1
fi
if [ ! -x "$(command -v gawk)" -o -z "$(command -v gawk)" ]; then
    echo -e "$helps\n!!! ERROR! Missing gawk.\n"
    exit 1
fi

# Read bedfile paths
shift $(($OPTIND - 1))
bedfiles=()
for bf in $*; do
    if [ -e $bf -a -n $bf ]; then
        bedfiles+=("$bf")
    else
        msg="!!!ERROR! Invalid bedfile, file not found.\n    File: $bf"
        echo -e " $helps\n$msg"
        exit 1
    fi
done
if [ 0 -eq ${#bedfiles[@]} ]; then
    msg="!!!ERROR! No bedfile was specified!\n"
    echo -e " $helps\n$msg"
    exit 1
fi

# Additional checks
if [ 1 -eq $csMode -a -z "$csList" -a 0 -eq $groupSize ]; then
    msg="!!!ERROR! Missing -l option with -c1 and -g."
    echo -e " $helps\n$msg"
    exit 1
fi
if [ ! $binStep -eq 0 -a $binSize -eq 0 ]; then
    echo -e "WARNING: missing bin size, ignoring -p option.\n"
fi
if [ ! $binSize -eq 0 -a $binStep -eq 0 ]; then
    binStep=$binSize
fi
if [ 0 -ne $binSize -a 0 -ne $groupSize ]; then
    if [ ! 0 -eq $(bc <<< "$binSize % $groupSize") ]; then
        msg="!!!ERROR! binSize ($binSize) must be divisible by groupSize."
        echo -e " $helps\n$msg"
        exit 1
    fi
fi
if [ 0 -ne $binStep -a 0 -ne $groupSize ]; then
    if [ $binStep -lt $groupSize ]; then
        msg="!!!WARNING! Using binStep as groupSize.\n            groupSize"
        msg="$msg should be smaller than or equal to binStep."
        echo -e " $helps\n$msg"
        groupSize=$binStep
    fi
fi

# Print settings
settings=""
if $chrWide; then
    settings="$settings
 Using chr-wide bins."
else
    settings="$settings
   Bin size : $binSize
   Bin step : $binStep"
fi

if [ 0 -ne $groupSize ]; then
    settings="$settings
 Group size : $groupSize"
fi

csModeLabel=()
csModeLabel+=("1:Universe")
csModeLabel+=("2:Union")
csModeLabel+=("3:Separate+Empty")
csModeLabel+=("4:Intersection")
settings="$settings
     Domain : ${csModeLabel[$csMode-1]}"
if [ 1 -eq $csMode ]; then
    settings="$settings
   Cutsites : $csList"
fi

if [ -n "$prefix" ]; then
    settings="$settings
     Prefix : '$prefix'"
fi

if [ -n "$suffix" ]; then
    settings="$settings
     Suffix : '$suffix'"
fi

if $normlast; then
    settings="$settings\n\n Normalizing over last condition."
fi

if $debugging; then
    settings="$settings\n\n Debugging mode ON."
fi

settings="$settings
 
 Output dir : $outdir
  Bed files :
   $(echo ${bedfiles[@]} | sed 's/ /\n   /g')"

echo -e "$settings\n"

# Set description
if $chrWide; then
    descr=$descr"bins.chrWide"
else
    descr=$descr"bins.size$binSize.step$binStep"
fi
if [ 0 -ne $groupSize ]; then
    descr="$descr.group$groupSize"
fi
descr="$descr.csm$csMode"
echo -e "$settings\n" > "$outdir/settings.$descr.txt"

# Constants
awkdir="`dirname ${BASH_SOURCE}`/awk/"

# RUN ==========================================================================

# 1) Identify chromosome sizes -------------------------------------------------
echo -e " Retrieving chromosome sizes ..."

# Retrieve chromosome sizes
chrSize=$(cat ${bedfiles[@]} | grep -v 'track' | datamash -sg1 -t$'\t' max 3)

# Sort by chromosomes
echo -e "$chrSize" | gawk -f "$awkdir/add_chr_id.awk" | sort -k1,1n | \
    cut -f2,3 > "$outdir/"$prefix"chr_size$suffix.tsv"


# 2) Generate bins -------------------------------------------------------------
echo -e " Generating bins ..."

# Generate bins
if $chrWide; then
    cat "$outdir/"$prefix"chr_size$suffix.tsv" | \
        gawk '{ print $1 "\t" 0 "\t" $2 }' \
        > "$outdir/"$prefix"$descr.bed" & pid=$!
else
    cat "$outdir/"$prefix"chr_size$suffix.tsv" | \
        gawk -v size=$binSize -v step=$binStep -f "$awkdir/mk_bins.awk" \
        > "$outdir/"$prefix"$descr.bed" & pid=$!
fi


# 3) Group reads ---------------------------------------------------------------

# Generate groups
if [ 0 -ne $groupSize ]; then
    echo -e " Generating groups ..."
    cat "$outdir/"$prefix"chr_size$suffix.tsv" | \
        gawk -v size=$groupSize -v step=$groupSize -f "$awkdir/mk_bins.awk" \
        > "$outdir/"$prefix"groups.$descr.bed" & pid=$!
fi

# Clean
wait $pid
if [ false == $debugging ]; then
    rm "$outdir/"$prefix"chr_size$suffix.tsv"
fi

if [ 0 -ne $groupSize ]; then
    echo -e " Grouping reads ..."

    # Group every bed file
    for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
        # Set input path
        infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{print $NF}')

        # Intersect with -loj to keep empty groups
        bedtools intersect -a "$outdir/"$prefix"groups.$descr.bed" \
            -b "${bedfiles[$bfi]}" -wa -wb -loj | cut -f 1-3,8 | \
            sed 's/-1$/0/' | gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
            > "$outdir/"$prefix"grouped.$descr.$infile$suffix.tsv" & pid=$!
        wait $pid

        # Point to group bed file instead of original one
        bedfiles[$bfi]="$outdir/"$prefix"grouped.$descr.$infile$suffix.tsv"
    done
fi


# 4) Normalize over last conditon ----------------------------------------------

if $normlast; then
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
fi
# Remove last condition
if [ false == $debugging ]; then rm "${bedfiles[-1]}"; fi
unset 'bedfiles[-1]'


# 5) Prepare and apply cutsite domain ------------------------------------------
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
            csbed=$(cat "$outdir/"$prefix"groups.$descr.bed" | grep -v "track")
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
            -a "$outdir/"$prefix"groups.$descr.bed" \
            -b <(echo -e "$csbed") -wa -c -loj | \
            cut -f 1-3 | sed 's/-1$/0/' | \
            gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
            > "$outdir/$prefix.cutsite.$descr.suffix.bed" & pid=$!
        wait $pid;
    else
        echo -e "$csbed" > "$outdir/$prefix.cutsite.$descr.suffix.bed"
    fi
fi

for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Set input/output paths (csd : CutSite Domain)
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
    if [ false == $normlast ]; then
        if [ 0 -ne $groupSize ]; then
            outfile=$(echo -e "$infile" | sed "s/grouped/csd/")
        else
            outfile=$prefix"csd.$descr$infile$suffix.tsv"
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

        bedtools intersect \
            -a "$outdir/"$prefix"groups.$descr.bed" \
            -b <(cat "$outdir/$outfile" | cut -f1-3) \
            -wa -c -loj | cut -f 1-3 | sed 's/-1$/0/' | \
            gawk -v prefix="row_" -f "$awkdir/add_name.awk" \
            > "$outdir/$prefix.cutsite.$descr.suffix.bed" & pid=$!
        wait $pid;
    fi

    # Remove grouped bed file
    if [ 0 -ne $groupSize -a false == $debugging ]; then
        rm "$outdir/$infile"
    fi

    # Point to non-zero-loci bed file instead of original one
    bedfiles[$bfi]="$outdir/$outfile"
done

# Clean
if [ false == $debugging ]; then
    rm "$outdir/"$prefix"groups.$descr.bed"
fi


# 6) Intersect with bedtools ---------------------------------------------------
echo -e " Assigning to bins ..."

# Assign bed reads to bins
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Input/output path
    infile=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
    outfile=$(echo -e "$infile" | sed 's/csd/intersected/')

    echo -e " > Assigning reads from $infile ..."
    bedtools intersect -a "$outdir/"$prefix"$descr.bed" \
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
    rm "$outdir/"$prefix"$descr.bed"
fi


# 7) Calculate bin statistics --------------------------------------------------
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
        
    fi

    # Remove binned
    wa/it $pid; if [ false == $debugging ]; then rm "${bedfiles[$bfi]}"; fi

    # Point to stats bed file instead of original one
    bedfiles[$bfi]="$outdir/$outfile"
done

exit 1
# 8) Assemble into bin data table ----------------------------------------------
echo -e " Combining information ..."

# combalize read count by cutsite and condition
comb=""
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
    # Combine
    cond_n_reads=$(cat "${bedfiles[$bfi]}" | grep -v "track" | datamash sum 4)
    tmp=$(cat "${bedfiles[$bfi]}" | \
        gawk -v cnr=$cond_n_reads -v bfi=$bfi -f "$awkdir/add_cnr_bfi.awk")
    comb="$comb$tmp\n"

    # Remove bin_stats
    wait $pid; if [ false == $debugging ]; then rm "${bedfiles[$bfi]}"; fi
done
# Columns of $comb:
# 1   2     3   4        5       6      7         8   9
# chr|start|end|condRead|readSum|readMu|readSigma|nCS|condID
comb=$(echo -e "$comb" | gawk -f "$awkdir/add_chr_id.awk" | \
    sort -k1,1n -k3,3n -k10,10n  | cut -f2- )

# Remove first line if empty
tmp=$(echo -e "$comb" | head -n1)
if [ -z "$tmp" ]; then comb=$(echo -e "$comb" | sed 1d); fi


# 9) Estimate centrality -------------------------------------------------------
echo -e " Estimating centrality ..."

# Prepare paste string
spaste=""; for i in $(seq 1 ${#bedfiles[@]}); do spaste="$spaste -"; done

#---------------------#
# Probability metrics #
#---------------------#

# Probability metric
echo -e " > Probability ..."
prob_mat=$(echo -e "$comb" | cut -f4,5,8 | \
    gawk -v type="p" -f "$awkdir/pre_process.awk" | paste $spaste)
probability_two_points=$(echo -e "$prob_mat" | \
    gawk -v calc="ratio" -v type="2p" -f "$awkdir/estimate_centrality.awk")
probability_fixed=$(echo -e "$prob_mat" | \
    gawk -v calc="ratio" -v type="f" -f "$awkdir/estimate_centrality.awk")
probability_global=$(echo -e "$prob_mat" | \
    gawk -v calc="ratio" -v type="g" -f "$awkdir/estimate_centrality.awk")

# Cumulative ratio metric
echo -e " > Cumulative ratio ..."
cumrat_mat="$prob_mat"
cumrat_two_points=$(echo -e "$cumrat_mat" | \
    gawk -v calc="ratio" -v cumrat=1 -v type="2p" \
    -f "$awkdir/estimate_centrality.awk")
cumrat_fixed=$(echo -e "$cumrat_mat" | \
    gawk -v calc="ratio" -v cumrat=1 -v type="f" \
    -f "$awkdir/estimate_centrality.awk")
cumrat_global=$(echo -e "$cumrat_mat" | \
    gawk -v calc="ratio" -v cumrat=1 -v type="g" \
    -f "$awkdir/estimate_centrality.awk")

# Ratio cumulative metric
echo -e " > Ratio cumulative ..."
ratcum_mat=$(echo -e "$comb" | cut -f4,5,8 | \
    tr '\t' ',' | paste $spaste)
ratcum_two_points=$(echo -e "$ratcum_mat" | \
    gawk -v calc="ratio" -v ratcum=1 -v type="2p" \
    -f "$awkdir/estimate_centrality.awk")
ratcum_fixed=$(echo -e "$ratcum_mat" | \
    gawk -v calc="ratio" -v ratcum=1 -v type="f" \
    -f "$awkdir/estimate_centrality.awk")
ratcum_global=$(echo -e "$ratcum_mat" | \
    gawk -v calc="ratio" -v ratcum=1 -v type="g" \
    -f "$awkdir/estimate_centrality.awk")

#---------------------#
# Variability metrics #
#---------------------#

# Variance metric
echo -e " > Variance ..."
var_mat=$(echo -e "$comb" | cut -f7 | paste $spaste)
var_two_points=$(echo -e "$var_mat" | \
    gawk -v calc="logratio" -v type="2p" -f "$awkdir/estimate_centrality.awk")
var_fixed=$(echo -e "$var_mat" | \
    gawk -v calc="logratio" -v type="f" -f "$awkdir/estimate_centrality.awk")

# Fano factor metric
echo -e " > Fano factor ..."
ff_mat=$(echo -e "$comb" | cut -f6,7 | \
    gawk -v type="ff" -f "$awkdir/pre_process.awk" | paste $spaste)
ff_two_points=$(echo -e "$ff_mat" | \
    gawk -v calc="diff" -v type="2p" -f "$awkdir/estimate_centrality.awk")
ff_fixed=$(echo -e "$ff_mat" | \
    gawk -v calc="diff" -v type="f" -f "$awkdir/estimate_centrality.awk")

# Coefficient of variation metric
echo -e " > Coefficient of variation ..."
cv_mat=$(echo -e "$comb" | cut -f6,7 | \
    gawk -v type="cv" -f "$awkdir/pre_process.awk" | paste $spaste)
cv_two_points=$(echo -e "$cv_mat" | \
    gawk -v calc="diff" -v type="2p" -f "$awkdir/estimate_centrality.awk")
cv_fixed=$(echo -e "$cv_mat" | \
    gawk -v calc="diff" -v type="f" -f "$awkdir/estimate_centrality.awk")

#----------#
# Assemble #
#----------#

# Prepare output table
metrics=$(echo -e "$comb" | cut -f1-3 | uniq | paste -d$'\t' - \
    <(echo -e "$probability_two_points") \
    <(echo -e "$probability_fixed") \
    <(echo -e "$probability_global") \
    <(echo -e "$cumrat_two_points") \
    <(echo -e "$cumrat_fixed") \
    <(echo -e "$cumrat_global") \
    <(echo -e "$ratcum_two_points") \
    <(echo -e "$ratcum_fixed") \
    <(echo -e "$ratcum_global") \
    <(echo -e "$var_two_points") \
    <(echo -e "$var_fixed") \
    <(echo -e "$ff_two_points") \
    <(echo -e "$ff_fixed") \
    <(echo -e "$cv_two_points") \
    <(echo -e "$cv_fixed") \
    )


# 10) Rank bins -----------------------------------------------------------------
echo -e " Ranking bins ..."

ranked=""
n_metrics=$(echo -e "$metrics" | awk '{print NF}' | sort -nu | tail -n 1)
for mi in $(seq 4 $n_metrics); do
    if $chrWide; then
        # Rank chromosomes and skip NaNs
        tmp=$(echo -e "$metrics" | cut -f1,$mi | gawk '"nan" != $4' | \
            sort -k2,2n | cut -f1)
    else
        # Rank sub-chromosome regions and skip NaNs
        tmp=$(echo -e "$metrics" | cut -f1,2,3,$mi | gawk '"nan" != $4' | \
            gawk -f "$awkdir/add_chr_id.awk" | \
            gawk '{ print $1"\t"$2"~"$3"~"$4"\t"$5 }' | \
            sort -k1,1n -k3,3n | cut -f2)
    fi
    if [ -z "$ranked" ]; then
        ranked=$tmp
    else
        ranked=$(paste -d$'\t' <(echo -e "$ranked") <(echo -e "$tmp"))
    fi
done


# 11) Output --------------------------------------------------------------------
echo -e " Writing output ..."

#------------#
# Add header #
#------------#

# combalization table
header="chr\tstart\tend"
header="$header\tcondRead\treadSum\treadMu\treadSigma\tnCS\tcondID"
comb=$(echo -e "$comb" | \
    gawk -v header="$header" -f "$awkdir/add_header.awk")

# Metrics table
header="chr\tstart\tend"
header="$header\tprob_2p\tprob_f\tprob_g"   # Probability
header="$header\tCoR_2p\tCoR_f\tCoR_g"      # Cumulative of Ratio
header="$header\tRoC_2p\tRoC_f\tRoC_g"      # Ratio of Cumulative
header="$header\tvar_2p\tvar_f"             # Variance
header="$header\tff_2p\tff_f"               # Fano factor
header="$header\tcv_2p\tcv_f"               # Coefficient of Variation
metrics=$(echo -e "$metrics" | \
    gawk -v header="$header" -f "$awkdir/add_header.awk")

# Ranked table
header=$(echo -e "$metrics" | head -n1 | cut -f4-)
ranked=$(echo -e "$ranked" | \
    gawk -v header="$header" -f "$awkdir/add_header.awk")

#-------#
# Write #
#-------#

# Remove bin positions if chromosome wide
if $chrWide; then
    comb=$(echo -e "$comb" | cut -f1,4-)
    metrics=$(echo -e "$metrics" | cut -f1,4-)
fi

# Write
echo -e "$comb" > "$outdir/"$prefix"combined.$descr$suffix.tsv"
echo -e "$metrics" > "$outdir/"$prefix"estimated.$descr$suffix.tsv"
echo -e "$ranked" > "$outdir/"$prefix"ranked.$descr$suffix.tsv"

# END ==========================================================================

################################################################################
