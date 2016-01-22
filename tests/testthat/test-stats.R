context("Exposure")

test_that("exposure calculations", {
  # Generating data
  set.seed(9999)
  diffnet <- rdiffnet(40,5, seed.p.adopt = .1)

  # Default
  exp_0_diffnet <- exposure(diffnet)
  exp_0_manual <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- diffnet$graph[[x]]
     ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(1e-15+Matrix::rowSums(s))
  })))

  expect_equivalent(exp_0_diffnet, exp_0_manual)

  # Struct Equiv
  exp_1_diffnet <- exposure(diffnet, wtype = 1)
  exp_1_manual <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- struct_equiv(diffnet$graph[[x]])$SE^(-1)
    s[!is.finite(s)] <- 0
    ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(1e-15+Matrix::rowSums(s))
  })))

  expect_equivalent(exp_1_diffnet, exp_1_manual)
})
