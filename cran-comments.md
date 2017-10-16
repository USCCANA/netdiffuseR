## Test environments

* Ubuntu 14.04 LTS (on travis-ci), R-rel 3.4.1 (2017-01-27), and R-rel 3.3.3 (2017-01-27).
* Windows Server 2012 R2 x64 (on AppVeyor), R version 3.4.2 (2017-09-28), R version 3.4.2 (2017-09-28).
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
