## Test environments

* Ubuntu 14.04.5 LTS 64-bit (on travis-ci), R 3.4.2, and R 3.3.3.
* Windows Server 2012 R2 (on AppVeyor), R version 3.4.2 (64-bit), R version 3.3.3 (32-bit).
* OS X El Capitan 10.11.6 (on travis-ci), R version 3.4.2, and R 3.3.3.

## R CMD check results

0 errors | 0 warnings | 1 note

* Checked with options --use-valgrind and --as-cran.
* Possible invalid URLs in man/rgraph_ws.Rd is false positive.
* Install size > 1Mb in libs. Using Rcpp increases the size. Still, working on
  reducing it.
* Also, possible mis-spelled words in DESCRIPTION have the correct spelling.

## Reverse dependencies

There are no reverse dependencies.
