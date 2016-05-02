## Test environments
* Local Ubuntu 14.04, R 3.2.5
* Ubuntu 12.04.5 LTS (on travis-ci), R version 3.2.5 (2016-04-14).
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
* About "Memtest notes: clang-UBSAN valgrind", the issue has been reported before
  (http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2014-December/008266.html)
  and it mostly arises from the usage of iterators in arma::sp_mat using
  RcppArmadillo.  I'm not aware of any solution for this (Valgrind pointing this
  issue), but the package includes several tests to make sure everything is in
  order (and will include more).
* Updated tests in testthat/*.R following the testthat package update (all passing).
  (Fixed errors that caused CRAN archiving the previous version).
  
