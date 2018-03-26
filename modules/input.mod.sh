#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 20170913
# Project: gpseq-centrality-estimate
# Description: this module manages help pages and input parameters.
# 
# ------------------------------------------------------------------------------



# MOD ==========================================================================

source "`dirname ${BASH_SOURCE}`/functions.mod.sh"

# Help string ------------------------------------------------------------------
helps="
 Usage: ./estimate_centrality.sh [-h][-d][-n][-c csMode][-l csList]
                                 [-b binList][-s binSize][-p binStep]
                                 [-g groupSize][-r prefix][-u suffix]
                                 [-e metrics][-i metrics]
                                 -o outdir [BEDFILE]...
"
long_helps="$helps
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

  # Select specific metrics
   By default, all the available metrics are calculated. Use -i to provide a
   list of the comma-separated metrics to calculate, while the rest would be
   excluded. Use the -e option to provide a list of the comma-separated metrics
   not to be calculated, while the rest would be included. The available metrics
   are:
    prob_2p prob_f prob_g
    cor_2p  cor_f  cor_g
    roc_2p  roc_f  roc_g
    var_2p  var_f
    ff_2p   ff_f
    cv_2p   cv_f

 Mandatory arguments:
  -o outdir     Output folder.
  BEDFILE       At least two (2) GPSeq condition bedfiles, space-separated and
                in increasing order of restriction conditions intensity.
                Expected to be ordered per condition. As BEDFILE is a positional
                argument, it should be provided after any other argument.

 Optional arguments:
  -h            Show this help page.
  -y            Do not ask for settings confirmation and proceed.
  -d            Debugging mode: save intermediate results.
  -n            Use last condition for normalization.
  -c csMode     Custite mode (see Notes).
  -l csList     Path to cutsite bed file. Required for -c1 when -g is not used.
  -b binList    path to bed file with bins. If used, -s and -p are ignored.
  -s binSize    Bin size in bp. Default to chromosome-wide bins.
  -p binStep    Bin step in bp. Default to bin sizeinStep.
  -g groupSize  Group size in bp. Used to group bins for statistics calculation.
                binSize must be divisible by groupSize. Not used by default.
  -e metrics    Comma-separated list of metrics to exclude from calculation.
                All metrics BUT the specified ones are calculated.
  -i metrics    Comma-separated list of metrics to be calculated.
                Only the specified metrics are calculated.
  -r prefix     Output name prefix.
  -u suffix     Output name suffix.
"

# Default values ---------------------------------------------------------------
csMode=3
binSize=0
binStep=0
groupSize=0
waitConfirm=true
chrWide=true            # If bins are chromosome wide
debugging=false         # If intermediate files are kept
normlast=false          # If normalization is performed
notOriginalBed=false    # If pointing to temporary bed file in $bedfiles

included_metrics=()
excluded_metrics=()
metrics="prob_2p,prob_f,prob_g,cor_2p,cor_f,cor_g,roc_2p,roc_f,roc_g,var_2p"
metrics="$metrics,var_f,ff_2p,ff_f,cv_2p,cv_f"
IFS=',' read -r -a metricslist <<< "$metrics"

# Parse options ----------------------------------------------------------------
while getopts hydnc:l:b:s:p:g:o:r:u:i:e: opt; do
    case $opt in
        h)
            # Help page
            echo -e "$long_helps" | less
            echo -e "$long_helps"
            exit 0
        ;;
        y)
            # Do not ask for confirmation
            waitConfirm=false
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
        b)
            # Bins bed file
            if [ -e $OPTARG -a -n "$OPTARG" ]; then
                binBed="$OPTARG"
                chrWide=false
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
        i)
            # Metrics to include
            IFS="," read -r -a included_metrics <<< "$OPTARG"
        ;;
        e)
            # Metrics to exclude
            IFS="," read -r -a excluded_metrics <<< "$OPTARG"
        ;;
        ?)
            msg="!!! ERROR! Unrecognized option."
            echo -e "$help\n$msg"
            exit 1
        ;;
    esac
done

# Check mandatory options ------------------------------------------------------
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

included=$(( 0 != ${#included_metrics[@]} ))
excluded=$(( 0 != ${#excluded_metrics[@]} ))
if (( 2 == $(bc <<< "$excluded + $included") )); then
  echo -e "$helps\n!!! ERROR! Options -e and -i cannot be used together."
  exit 1
fi
if (( 1 == $included )); then
  calc_metrics=()
  for metlab in ${included_metrics[@]}; do
    containsElement "$metlab" "${metricslist[@]}"
    if [[ 0 == "$?" ]]; then
      echo -e "$helps\n!!! ERROR! Metric '$metlab' not recognized (-i)."
      echo -e "Allowed metrics: ${metricslist[@]}"
      exit 1
    fi
  done
  for metlab in ${metricslist[@]}; do
    containsElement "$metlab" "${included_metrics[@]}"
    if [[ 1 == "$?" ]]; then
      calc_metrics+=("$metlab")
    fi
  done
fi
if (( 1 == $excluded )); then
  calc_metrics=${metricslist[@]}
  for metlab in ${excluded_metrics[@]}; do
    containsElement "$metlab" "${metricslist[@]}"
    if [[ 0 == "$?" ]]; then
      echo -e "$helps\n!!! ERROR! Metric '$metlab' not recognized (-e)."
      echo -e "Allowed metrics: ${metricslist[@]}"
      exit 1
    else
      calc_metrics=("${calc_metrics[@]/$metlab}")
      calc_metrics=$(join_by "," ${calc_metrics[@]})
      calc_metrics=($(split_by "," "$calc_metrics"))
    fi
  done
fi
if (( 0 == $(bc <<< "$excluded + $included") )); then
  calc_metrics=${metricslist[@]}
fi

# Read bedfile paths -----------------------------------------------------------
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

# Additional checks ------------------------------------------------------------
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

# Prepare to print settings ----------------------------------------------------
settings=" # GPSeq-centrality-estimate\n"
if [ -n "$binBed" ]; then
  settings="$settings\n Bin bed file: $binBed"
else
  if $chrWide; then settings="$settings\n Using chr-wide bins.";
  else settings="$settings\n   Bin size : $binSize\n   Bin step : $binStep"; fi
fi

if [ 0 -ne $groupSize ]; then settings="$settings\n Group size : $groupSize"; fi

csModeLabel=()
csModeLabel+=("1:Universe")
csModeLabel+=("2:Union")
csModeLabel+=("3:Separate+Empty")
csModeLabel+=("4:Intersection")
settings="$settings\n     Domain : ${csModeLabel[$csMode-1]}"
if [ 1 -eq $csMode ]; then settings="$settings\n   Cutsites : $csList"; fi

if [ -n "$prefix" ]; then settings="$settings\n     Prefix : '$prefix'"; fi
if [ -n "$suffix" ]; then settings="$settings\n     Suffix : '$suffix'"; fi
if $normlast; then settings="$settings\n\n Normalizing over last condition."; fi
if $debugging; then settings="$settings\n\n Debugging mode ON."; fi

if(( 1 == $included )); then 
  settings="$settings\n\n Included metrics:\n"
  settings="$settings   $(join_by ", " ${included_metrics[@]} | sed 's/,/, /g')"
fi
if(( 1 == $excluded )); then 
  settings="$settings\n\n Excluded metrics:\n"
  settings="$settings   $(join_by ", " ${excluded_metrics[@]} | sed 's/,/, /g')"
fi

settings="$settings\n 
 Output dir : $outdir
  Bed files :"
for bfi in $(seq 0 ${#bedfiles[@]}); do
  if [ -n "${bedfiles[bfi]}" ]; then
    bfip=$(bc <<< "$bfi+1")
    settings="$settings\n   $bfip : ${bedfiles[bfi]}"
  fi
done

# Ask for confirmation ---------------------------------------------------------
end=1

if $waitConfirm; then
  settings_confirm="
   ##############################################
   #                                            #
   #  PLEASE, DOUBLE CHECK THE SETTINGS BELOW   #
   # (press 'q' to continue when you are ready) #
   #                                            #
   ##############################################

  $settings\n\n"
  echo -e "$settings_confirm" | less
  msg="$msg\nRun the analysis?\nYes (y), Abort (a), Show again (s)"
  clear; echo -e $msg; read -e ans
  end=0
fi
while [[ 0 -eq $end ]]; do
  if [[ -z $ans ]]; then
    echo -e $msg
    read -e ans
  elif [[ 'a' == $ans ]]; then
    end=1
    echo "Aborted."
    exit 1
  elif [[ 'y' == $ans ]]; then
    echo -e "\n"
    end=1
  elif [[ 's' == $ans ]]; then
    echo -e "$settings_confirm" | less
    clear; echo -e $msg; read -e ans
  else
    echo -e $msg
    read -e ans
  fi
done

# Re-print to display before proceeding
clear; echo -e "\n$settings\n"

# Save settings to file --------------------------------------------------------
if [ -n "$binBed" ]; then
  descr=$descr"customBins"
else
  if $chrWide; then descr=$descr"bins.chrWide";
  else descr=$descr"bins.size$binSize.step$binStep"; fi
fi
if [ 0 -ne $groupSize ]; then descr="$descr.group$groupSize"; fi
descr="$descr.csm$csMode"
if $normlast; then descr="$descr.norm"; fi
echo -e "$settings\n" > "$outdir/$prefix""settings.$descr.txt"

# END ==========================================================================

################################################################################
