## Test environments

* Local Ubuntu 14.04, R-rel 3.3.2 (2016-10-31).
* Ubuntu 12.04.5 LTS (on travis-ci), R-rel 3.3.2 (2016-06-21).
* Windows Server 2012 R2 x64 (on AppVeyor), R-dev 2016-11-01 r71616, and R-rel 3.3.2 (2016-10-31).
* OS X El Capitan 10.11.6 (on travis-ci), R-rel version 3.3.2 (2016-10-31).

## R CMD check results

0 errors | 0 warnings | 1 note

* Checked with options --use-valgrind and --as-cran.
* Possible invalid URLs in man/rgraph_ws.Rd is false positive.
* Install size > 1Mb in libs. Using Rcpp increases the size. Still, working on
  reducing it.
* Also, possible mis-spelled words in DESCRIPTION have the correct spelling.

## Reverse dependencies

There are no reverse dependencies.
