# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).



## Unreleased
### ...
- ...



## [2.3.5] - 2019-03-25
### Fixed
- `gpseqc_estimate`
    + Fixed `-l` option.



## [2.3.4] - 2018-09-05
### Fixed
- `gpseqc_estimate`
    + Settings confirmation page.
    + Removed NaN when calculating IQR and quartiles. Discarding NaN when asking for outliers with IQR method and `lim` specified.
    + Colliding temporary files issue using `tempfile.TemporaryDirectory()`.



## [2.3.3] - 2018-08-16
### Added
- `gpseqc_estimate`
    + Additional help page with `-H`.
    + Added "rescaled" scores output, with scores rescaled based as `2**((score-min(score)/max(score-min(score))))`. Only non-outlier scores are considered in the calculation of min/max. Hence, non-outlier scores are in the [1;2] interval, lower outliers in the [0;1) interval, and upper outliers in the (2;Inf] interval.
    + Added parameters to define outliers for rescaling: `--score-outliers`, `--score-outlier-alpha` and `--score-outlier-lim`.

### Changed
- `gpseqc_estimate`
    + Use only chromosomes present in the input bed files when running chromosome-wide with `-G` option.
    + Split script help page in `-h` for attributes and standard help, and `-H` for more descriptive and readable text.
    + Renamed bed outlier parameters to avoid confusion with new score outliers:
        * `--outliers` to `--bed-outliers`
        * `--outlier-alpha` to `--bed-outlier-alpha`
        * `--outlier-lim` to `--bed-outlier-lim`
    + Renamed second parameter of `gpseqc.stats.score_outliers` to `method` (previously `stype`).



## [2.3.2] - 2018-08-16
## Fixed
- `gpseqc_estimate`
    + Fixed masking with `-M`, now removing all bins overlapping of even 1 bp with the masked regions.



## [2.3.1] - 2018-07-09
### Added
- Option to mask `gpseqc_estimate` output (`-M`).



## [2.3.0] - 2018-06-26
### Added
- `gpseqc_estimate`
    + New variability-based estimates, where conditions are combined in the same way as the probability-based ones. The names are the same, with a prefix `s`.
    + Estimates list and ermetic description in main help page.



## [2.2.0] - 2018-06-07
### Added
- `gpseqc_estimate`
    + Outlier export (only in `debug` mode).
    + Added `alpha/lim` outlier options.
    + Forced flushing of temporary files.
    + Fixed *sorting* issues.
    + Bed file masking at the beginning of the pipeline.

### Changed
- `gpseqc_estimate`
    + Remove all outliers as default behavior.



## [2.1.0] - 2018-06-04
### Added
- Outlier filter to `gpseqc_estimate`.



## [2.0.7] - 2018-05-26
### Fixed
- Axis drop bug in `gpseqc_compare`.
- Matplotlib back-end to work when DISPLAY is unavailable.
- `-1` bug when intersecting in `gpseqc_estimate`.



## [2.0.6] - 2018-05-22
###  Fixed
- Assert error triggered by Kendall tau weighted when all elements in the ranking have the same weight.

### Added
- Option to allow for custom centrality estimates when comparing rank tables with `gpseqc_compare` (`-C/--custom-estimates`).
- Chromosome size option to `gpseqc_estimate` when using chromosome wide bins. Using genome size file as in `bedtools`.



## [2.0.5] - 2018-05-21
### Fixed
- Sorting bug in `gpseqc_compare`.
- EMD calculation when `np.nan` values are present, in `gpseqc_compare`.



## [2.0.4] - 2018-05-16
### Added
- `gpseqc_compare` for ranking comparison.



### [2.0.3] - 2018-05-09
### Added
- `CMETRICS` constant to `gpseqc_estimate` with lambda functions to calculate centrality metrics.
- `RankTable` class to `comparison.py`.
- Added timestamp, command line and version tag to `gpseqc_estimate` settings file for easier bug reporting and debugging.

### Fixed
- `-y` option now working in `gpseqc_estimate`.
- Typo in `gpseqc_estimate`.



## [2.0.2] - 2018-04-24
### Fixed
- Empty bin issue during step 6 of `gpseqc_estimate`.
- Division by zero in `centrality.py`.



## [2.0.1] - 2018-04-24
### Added
- Error triggered when empty bed files are provided to `gpseqc_estimate`.



## [2.0.0] - 2018-04-24
### Changed
- Implemented as Python package `gpseqc`.
- Rank comparison still available only as legacy version.
- Documentation moved to [GitHub Pages](https://ggirelli.github.io/gpseq_ce/).



## [1.2.2] - 2018-04-16
### Fixed
- Settings file output name now includes suffix.



## [1.2.1] - 2018-03-27
### Fixed
- Fixed behavior when -i/-e are not used.

### Added
- Version option.



## [1.2.0] - 2018-03-23
### Added
- Output directory parameter to comparison script.
- Option to include/exclude metrics from calculation.

### Fixed
- Comparison script displays error message when bed file is required but not provided.
- Now allowing to plot when DISPLAY is disconnected.

### Changed
- Fixed order of settings in confirmation message.
- Added common functions in separate module.



## [1.1.0] - 2018-02-13
### Added
- Implemented rank comparison script.

### Fixed
- Not skipping first bed file in settings confirmation page.
- Added prefix to output settings file.



## [1.0.0] - 2018-01-17 - First release.



[Unreleased] https://github.com/ggirelli/gpseq-centrality-estimate  
[2.3.4] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.3.4  
[2.3.3] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.3.3  
[2.3.2] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.3.2  
[2.3.1] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.3.1  
[2.3.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.3.0  
[2.2.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.2.0  
[2.1.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.1.0  
[2.0.7] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.7  
[2.0.6] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.6  
[2.0.5] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.5  
[2.0.4] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.4  
[2.0.3] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.3  
[2.0.2] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.2  
[2.0.1] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.1  
[2.0.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.0  
[1.2.2] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.2.2  
[1.2.1] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.2.1  
[1.2.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.2.0  
[1.1.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.1.0  
[1.0.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.0.0  