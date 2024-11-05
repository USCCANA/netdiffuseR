
test_that("multidiffusion exposure calculations", {
  # Generating data
  diffnet <- rdiffnet(40,5, seed.p.adopt = .1)

  #data(medInnovationsDiffNet)
  #exposure(medInnovationsDiffNet)

  # two spreads
  cumadopt_2 <- medInnovationsDiffNet$cumadopt
  cumadopt_2 <- array(c(cumadopt_2,cumadopt_2[rev(1:nrow(cumadopt_2)),]), dim=c(dim(cumadopt_2), 2))

  # Default
  #ans0 <- exposure(medInnovationsDiffNet)#exposure(diffnet)
  ans1 <- as.matrix(do.call(cbind,lapply(medInnovationsDiffNet$meta$pers, function(x) {
    s <- medInnovationsDiffNet$graph[[x]]
    for (q in 1:dim(cumadopt)[3]) {
      ( s %*% cumadopt_2[,x,q,drop=FALSE])/(1e-15+Matrix::rowSums(s))
    }
  })))

  exposure(medInnovationsDiffNet$graph, medInnovationsDiffNet$cumadopt)

  ans2 <- exposure(medInnovationsDiffNet$graph, cumadopt = cumadopt_2)
  #ans3 <- exposure(as.array(medInnovationsDiffNet), cumadopt = cumadopt_2)

  expect_equivalent(ans0, ans1)
  expect_equivalent(ans0, ans2)
  expect_equivalent(ans0, ans3)

  # With an attribute
  X <- matrix(diffnet[["real_threshold"]], ncol=5, nrow=40, byrow = FALSE)
  ans0 <- exposure(diffnet, attrs=X)
  ans1 <- exposure(diffnet, attrs="real_threshold")
  expect_equivalent(ans0, ans1)

  expect_error(exposure(diffnet$graph, attrs="real_threshold"),"is only valid for")

  # Struct Equiv
  se <- struct_equiv(diffnet)
  se <- lapply(se, function(x) {
    ans <- methods::as(x$SE, "dgCMatrix")
    ans@x <- 1/(ans@x + 1e-20)
    ans
  })
  exp_1_diffnet <- exposure(diffnet, alt.graph = se, valued=TRUE)
  se2 <- vector("list", length(se))
  exp_1_manual <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- methods::as(struct_equiv(diffnet$graph[[x]])$SE, "dgCMatrix")
    s@x <- 1/(s@x + 1e-20)
    se2[[x]] <<- s
    ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(Matrix::rowSums(s) +1e-20)
  })))

  expect_equivalent(unname(exp_1_diffnet), unname(exp_1_manual))

  # Lagged exposure
  ans0 <- exposure(diffnet)
  ans1 <- exposure(diffnet, lags = 1)
  ans2 <- exposure(diffnet, lags = 2)
  ans3 <- exposure(diffnet, lags = -1)

  expect_equivalent(ans0[,-5], ans1[,-1])
  expect_equivalent(ans0[,-(4:5)], ans2[,-(1:2)])
  expect_equivalent(ans0[,-1], ans3[,-5])

  expect_error(exposure(diffnet, lags=5), "cannot be greater")
  expect_error(exposure(diffnet, lags=NA))
  expect_error(exposure(diffnet, lags=c(1:2)))

})

test_that("multidiffusion exposure calculations", {
  # Generating data
  set.seed(999)
  diffnet <- rdiffnet(40,5, seed.p.adopt = .1)

  # Default
  ans0 <- exposure(diffnet)
  ans1 <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- diffnet$graph[[x]]
    ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(1e-15+Matrix::rowSums(s))
  })))
  ans2 <- exposure(diffnet$graph, cumadopt = diffnet$cumadopt)
  ans3 <- exposure(as.array(diffnet), cumadopt = diffnet$cumadopt)

  expect_equivalent(ans0, ans1)
  expect_equivalent(ans0, ans2)
  expect_equivalent(ans0, ans3)

  # With an attribute
  X <- matrix(diffnet[["real_threshold"]], ncol=5, nrow=40, byrow = FALSE)
  ans0 <- exposure(diffnet, attrs=X)
  ans1 <- exposure(diffnet, attrs="real_threshold")
  expect_equivalent(ans0, ans1)

  expect_error(exposure(diffnet$graph, attrs="real_threshold"),"is only valid for")

  # Struct Equiv
  se <- struct_equiv(diffnet)
  se <- lapply(se, function(x) {
    ans <- methods::as(x$SE, "dgCMatrix")
    ans@x <- 1/(ans@x + 1e-20)
    ans
  })
  exp_1_diffnet <- exposure(diffnet, alt.graph = se, valued=TRUE)
  se2 <- vector("list", length(se))
  exp_1_manual <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- methods::as(struct_equiv(diffnet$graph[[x]])$SE, "dgCMatrix")
    s@x <- 1/(s@x + 1e-20)
    se2[[x]] <<- s
    ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(Matrix::rowSums(s) +1e-20)
  })))

  expect_equivalent(unname(exp_1_diffnet), unname(exp_1_manual))

  # Lagged exposure
  ans0 <- exposure(diffnet)
  ans1 <- exposure(diffnet, lags = 1)
  ans2 <- exposure(diffnet, lags = 2)
  ans3 <- exposure(diffnet, lags = -1)

  expect_equivalent(ans0[,-5], ans1[,-1])
  expect_equivalent(ans0[,-(4:5)], ans2[,-(1:2)])
  expect_equivalent(ans0[,-1], ans3[,-5])

  expect_error(exposure(diffnet, lags=5), "cannot be greater")
  expect_error(exposure(diffnet, lags=NA))
  expect_error(exposure(diffnet, lags=c(1:2)))

})
