Estimate centrality
===

The **estimate centrality** tool aims at estimating the 3D spatial nuclear centrality of genomic regions.

### estimate_centrality.sh

Two main families of metrics are implemented: probability-based and variance-based.

Added empty loci removal and normalization over last condition. Fixed intermediate file naming and output. Counted cutsites are now actually seen custites, disregarding empty cutsites (with no reads).

```
usage: ./estimate_centrality.sh [-h][-d][-s binSize][-p binStep][-g groupSize]
                                [-r prefix][-u suffix] -o outdir [BEDFILE]...

 Description:
  Estimate global centrality. The script performs the following steps:
   (1) Identify & sort chromosomes
   (2) Generate bins
   (3) Group cutsites (intersect)
   (4) Remove empty cutsites/groups
   (5) Normalize over last condition.
   (6) Assign grouped reads to bins (intersect)
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
  Statistics (mean, variance) metrics take into account only cutsites sensed
  in that condition. The script ignores 'zero' loci (with no reads). This is
  true for both probability- and variability-based metrics.

  Depending on the sequencing resolution, it might not be feasible to go for
  single-cutsite resolution. Thus, cutsite can be grouped for the statistics
  calculation using the -g option.

  In case of sub-chromosome bins, the ranking is done in an ordered
  chromosome-wise manner.

  Empty custites/groups are removed before bin assignment, while empty bins are
  kept. Also, normalization is performed after empty bin removal but before bin
  assignment, i.e., either on the grouped or single cutsites.

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
  -s binSize    Bin size in bp. Default to chromosome-wide bins.
  -p binStep    Bin step in bp. Default to bin sizeinStep.
  -g groupSize  Group size in bp. Used to group bins for statistics calculation.
                binSize must be divisible by groupSize. Not used by default.
  -r prefix     Output name prefix.
  -u suffix     Output name suffix.
```

### rankcmp.py

A script to compare rankings of the same type. It is compatible with both FISH and Sequencing based rankings. The rank type can be one of `chr`, `set`, `probe` for compatibility with the FISH-based input. When comparing two seq-based ranks, use `chr` for chrWide binning and either `set` or `probe` for sub-chromosome bins.

The following Python2 packages should be installed beforehand: `argparse`, `joblib`, `math`, `matplotlib`, `matplotlib`, `numpy`, `pandas`, `tqdm`, `random`, `scipy`, `time`, `warnings`.

```
usage: rankcmp.py [-h] [-t type] [-d distance] [--bed path] [-i niter]
                  [--threads nthreads] [-s delimiter] [--prefix text]
                  rank1 rank2

Calculate difference between rankings obtained with GPSeq, either by FISH or
sequencing.

positional arguments:
  rank1               Path to first ranking set.
  rank2               Path to second ranking set.

optional arguments:
  -h, --help          show this help message and exit
  -t type             One of the following: chr, set, probe.
  -d distance         Wither 'kt' (Kendall tau) or 'ktw' (Kendall tau
                      weighted).
  --bed path          Path to bed file.
  -i niter            Number of iterations to build the random distribution.
  --threads nthreads  Number of threads for parallelization.
  -s delimiter        Ranking file field delimiter.
  --prefix text       Text for output file name prefix.
```