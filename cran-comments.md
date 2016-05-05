## Test environments
* Local Ubuntu 14.04, R version 3.3.0 (2016-05-03)
* Ubuntu 12.04.5 LTS (on travis-ci), R version 3.3.0 (2016-05-03).
* Windows Server 2012 R2 x64 (on AppVeyor), R Under development (3.2.5 Patched (2016-04-18 r70556).
* OS X 10.9.5 (Mavericks) (on travis-ci), R version 3.2.4 (2016-03-10).

## R CMD check results

0 errors | 0 warnings | 1 note

* Install size > 1Mb in libs. Using Rcpp increases the size. Still, working on
  reducing it.
* Also, possible mis-spelled words in DESCRIPTION have the correct spelling.
* CRAN repository db overrides: X-CRAN-Comment: Archived on 2016-05-02 as ...
  fixed (see below).
* Possible invalid URLs in man/rgraph_ws.Rd is false positive.

## Reverse dependencies

There are no reverse dependencies.

## Resubmission

This is a resubmission. In this version I have:

* Included rgb (grDevices), polygon (graphics) and setNames (stats) as importFrom
  in the package's NAMESPACE.
* Fixed the problem with "Vignette sources missing".
* About "Memtest notes: clang-UBSAN valgrind", the issue was actually a bug in
  RcppArmadillo. This has been solved in RcppArmadillo v 0.6.700.6.0 (which is
  already on CRAN). R CMD check --use-valgrind no longer reports memory leaks
  or 'Conditional jump...' issues.
* Updated tests in testthat/*.R following the testthat package update (all passing).
  (Fixed errors that caused CRAN archiving the previous version).
  
