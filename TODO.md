
# gpseqc_estimate

* Fix normalization to take number of cutsites into consideration.
* Manually verify that cs_mode != 3 works as expected.
* When running in chomosome-wide mode with `-G`, use only chromosomes present in the input bed files.
* Implement score normalization for inter-sample comparison.

# gpseqc_compare

* Move `nan` removal outside of `compare()` to save time during iterations.
