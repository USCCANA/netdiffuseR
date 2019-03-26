* This is a resubmission: Fixing broken link in manual.

## R CMD check results

* Checked with options --use-valgrind and --as-cran.

* Possible invalid URLs in man/rgraph_ws.Rd is false positive.

* Install size > 1Mb in libs. Using Rcpp increases the size. Still, working on
  reducing it.

* Also, possible mis-spelled words in DESCRIPTION have the correct spelling.

* Error in R Old Rel (3.4.4) has to do with dependency change in statnet.common
  # which I cannot control.
  
* One known error while installing in OldRel because of a change in the system
  requirements of the package statnet.common (on which netdiffuseR depends).

## Reverse dependencies

There are no reverse dependencies.
