# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).



## Unreleased
### ...
- ...



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
[2.0.2] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.2  
[2.0.1] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.1  
[2.0.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v2.0.0  
[1.2.2] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.2.2  
[1.2.1] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.2.1  
[1.2.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.2.0  
[1.1.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.1.0  
[1.0.0] https://github.com/ggirelli/gpseq-centrality-estimate/releases/tag/v1.0.0  