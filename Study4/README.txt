The main analyses are run by running: Pgm/Ana_02DEC2014.R.This script calls other scripts:
	- new_function.R includes the custom function needed
	- code_variables_cat.R QC the data and categorise continuous variables
	- impute_uk_bio*.R impute data using mice package
	- prediction*.R run lasso regression
	- univariate_noint*.R univariate association analysis NON time-dependent
	- univariate_yesint1*.R univariate association analysis time-dependent, only association results
	- univariate_yesint2*.R univariate association analysis time-dependent, only C-index results

Several of these scripts can be run in parallel using a computer cluster