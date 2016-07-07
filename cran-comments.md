## Test environments

* Local Ubuntu 14.04, R version 3.3.1 (2016-06-21).
* Ubuntu 12.04.5 LTS (on travis-ci), version 3.3.1 (2016-06-21).
* Windows Server 2012 R2 x64 (on AppVeyor), R version 3.3.1 Patched (2016-06-28 r70858).
* OS X 10.9.5 (Mavericks) (on travis-ci), R version 3.3.1 (2016-06-21).

## R CMD check results

0 errors | 0 warnings | 1 note

* Possible invalid URLs in man/rgraph_ws.Rd is false positive.
* Install size > 1Mb in libs. Using Rcpp increases the size. Still, working on
  reducing it.
* Also, possible mis-spelled words in DESCRIPTION have the correct spelling.

## Reverse dependencies

There are no reverse dependencies.

## Resubmission

This is a resubmission. In this version I have:

* Reduced the times on examples, tests and vignettes as requested by CRAN.
  
