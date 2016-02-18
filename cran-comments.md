## Test environments
* local Ubuntu 14.04 install, R 3.2.3
* Ubuntu 12.04.5 LTS (on travis-ci), R version 3.2.3 (2015-12-10)
* Windows Server 2012 R2 x64 (on AppVeyor), R Under development (unstable) (2016-02-15 r70179)

## R CMD check results

0 errors | 0 warnings | 2 note

* This is a new release.
* Install size > 1Mb in libs. Using Rcpp increases the size. Still, working on
  reducing it.

* Also, possible mis-spelled words in DESCRIPTION have the correct spelling.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

## Resubmission
This is a resubmission. In this version I have:

* Included DOI/ISBN references in DESCRIPTION
* Fixed the warnings:
  ' infection.cpp:35:20: warning: ISO C++ forbids variable length array ‘discount’ [-Wvla]'
  and
  ' infection.cpp:136:20: warning: ISO C++ forbids variable length array ‘discount ’ [-Wvla]'

