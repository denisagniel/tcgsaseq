# News about the `tcgsaseq` R package


### Main changes in Version 1.8.0 (2017-11-23):
* dsFDR function for accurate discrete FDR control when using permutation test. This 
changes the column names of the `$pvals` output when `which_test = permutations`


### Main changes in Version 1.7.2 (2017-07-18):
* faster permutations

### Main changes in Version 1.7.1 (2017-07-07):
* faster implementation

### Main changes in Version 1.6.5 (2017-12-24) --- *this is only a minor release*:
* bug fix in the log2-cpm transformation computation (for `preprocessed = FALSE`) in `sp_weights`
* bug fix in `sp_weights` when facing `NA` or `NaN`values


### Main changes in Version 1.6.3 (2017-11-09)
* `WARNING` instead of `ERROR` when a gene set with no measured genes is tested
* `NA` support through `na.rm_...` logical arguments
* bug fix for row.names error in weights computation for gene-wise analysis


### Main changes in Version 1.6.2 (2017-08-04)
* bug fix to the log2-cpm transformation computation (for `preprocessed = FALSE`)
* improvement to defaults in `varseq` and `tcgsa_seq` wrapper functions, including
 an option to compute heteroscedasticity weights without conditioning on the variable(s)
 to be tested


### Main changes in Version 1.5.2 (2017-07-11)
* user-friendly wrapper function `varseq` for gene-wise testing
* RAM usage optimization


### Main changes in Version 1.5.1 (2017-05-24) --- *this is only a minor release*:
* homogeneous test now available for gene-wise testing.


### Main changes in Version 1.5.0 (2017-03-27):
* permutation test now available for gene-wise testing


### Main changes in Version 1.4.0 (2017-02-21):
* bug fix for gene-wise testing


### Main changes in Version 1.3.1 (2016-11-24):
* bug fix for testing several covariates together (such as several bases of time) which now gives correct p-values
* support for gene-wise testing
* example dataset 'baduel_small' was changed to 'baduel_5gs' which includes 1943 additional gene expressions and the definition of 5 gene sets
* documentation updated


### Main changes in Version 1.2.0 (2016-07-19):
* support homogeneous gene set test
* documentation updated


### Main changes in Version 1.1.0 (2016-05-26):
* improved imports
* documentation updated


### Main changes in Version 1.0.2 (2016-05-18) --- *this is only a minor release*:
* documentation and CITATION updated


### Main changes in Version 1.0.1 (2016-05-14) --- *this is only a minor release*:
* help file improvements

