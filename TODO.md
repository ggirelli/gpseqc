
# gpseqc_estimate

* Fix normalization to take number of cutsites into consideration.
* Manually verify that cs_mode != 3 works as expected.

# gpseqc_compare

* Move `nan` removal outside of `compare()` to save time during iterations.
