library(testthat)
library(netdiffuseR)

Sys.setenv("R_TESTS" = "")

test_check("netdiffuseR", reporter = "summary")
