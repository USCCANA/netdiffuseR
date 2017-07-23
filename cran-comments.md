## Resubmission

* While version 1.18.0 was recently sent to CRAN, this version fixes
  some critical bugs that were identified right after that.
  

## Test environments

* Local Ubuntu 14.04.5 LTS, R version 3.4.0 (2017-04-21).
* Ubuntu 12.04.5 LTS (on travis-ci), R-rel 3.4.1 (2017-06-30).
* Windows Server 2012 R2 x64 (on AppVeyor), R Under development (unstable) (2017-07-17 r72924), and R version 3.4.1 (2017-06-30).
* OS X El Capitan 10.11.6 (on travis-ci), R version 3.4.1 (2017-06-30).

## R CMD check results

0 errors | 0 warnings | 1 note

* Checked with options --use-valgrind and --as-cran.
* Possible invalid URLs in man/rgraph_ws.Rd is false positive.
* Install size > 1Mb in libs. Using Rcpp increases the size. Still, working on
  reducing it.
* Also, possible mis-spelled words in DESCRIPTION have the correct spelling.

## Reverse dependencies

There are no reverse dependencies.
