test_that("multidiffusion exposure calculations", {
  # Generating data
  diffnet <- rdiffnet(40,5, seed.p.adopt = .1)

  # Creating two spreads
  cumadopt_2 <- diffnet$cumadopt
  cumadopt_2 <- array(c(cumadopt_2,cumadopt_2[rev(1:nrow(cumadopt_2)),]), dim=c(dim(cumadopt_2), 2))

  # Default --
  ans0 <- exposure(diffnet, cumadopt = cumadopt_2)
  ans1 <- array(unlist(lapply(1:dim(cumadopt_2)[3], function(q) {
    lapply(diffnet$meta$pers, function(x) {
      graph_slice <- diffnet$graph[[x]]
      as.numeric((graph_slice %*% cumadopt_2[, x, q, drop = FALSE]) /
                   (1e-20 + Matrix::rowSums(graph_slice)))
    })
  })), dim = dim(cumadopt_2))

  ans2 <- exposure(diffnet$graph, cumadopt = cumadopt_2)
  ans3 <- exposure(as.array(diffnet), cumadopt = cumadopt_2)

  #ans0 - ans1
  expect_equivalent(ans0, ans1)
  expect_equivalent(ans0, ans2)
  expect_equivalent(ans0, ans3)

  # By each behavior --
  ans4 <- exposure(diffnet)
  ans5 <- exposure(diffnet$graph, cumadopt = diffnet$cumadopt)
  cumadopt_rev <- diffnet$cumadopt[rev(1:nrow(diffnet$cumadopt)),]
  ans6 <- exposure(diffnet$graph, cumadopt = cumadopt_rev)

  expect_equivalent(ans0[,,1], ans4)
  expect_equivalent(ans0[,,1], ans5)
  expect_equivalent(ans0[,,2], ans6)

  # With an attribute --
  X <- matrix(diffnet[["real_threshold"]], ncol=5, nrow=40, byrow = FALSE)
  ans0 <- exposure(diffnet$graph, cumadopt = cumadopt_2, attrs=X)
  ans1 <- exposure(as.array(diffnet), cumadopt = cumadopt_2, attrs=X)
  expect_equivalent(ans0, ans1)

  expect_error(exposure(diffnet$graph, attrs="real_threshold"),"is only valid for")

  # Struct Equiv --
  se <- struct_equiv(diffnet)
  se <- lapply(se, function(x) {
    ans <- methods::as(x$SE, "dgCMatrix")
    ans@x <- 1/(ans@x + 1e-20)
    ans
  })
  ans0 <- exposure(diffnet, cumadopt = cumadopt_2, alt.graph = se, valued=TRUE)
  ans1 <- array(unlist(lapply(1:dim(cumadopt_2)[3], function(q) {
    lapply(diffnet$meta$pers, function(x) {
      graph_slice <- methods::as(struct_equiv(diffnet$graph[[x]])$SE, "dgCMatrix")
      graph_slice@x <- 1/(graph_slice@x + 1e-20)
      as.numeric((graph_slice %*% cumadopt_2[, x, q, drop = FALSE]) /
                   (1e-20 + Matrix::rowSums(graph_slice)))
    })
  })), dim = dim(cumadopt_2))

  #ans0 - ans1
  expect_equivalent(unname(ans0), unname(ans1))

  # Lagged exposure --
  ans0 <- exposure(diffnet, cumadopt = cumadopt_2)
  ans1 <- exposure(diffnet, cumadopt = cumadopt_2, lags = 1)
  ans2 <- exposure(diffnet, cumadopt = cumadopt_2, lags = 2)
  ans3 <- exposure(diffnet, cumadopt = cumadopt_2, lags = -1)

  expect_equivalent(ans0[,-5,], ans1[,-1,])
  expect_equivalent(ans0[,-(4:5),], ans2[,-(1:2),])
  expect_equivalent(ans0[,-1,], ans3[,-5,])

  expect_error(exposure(diffnet, lags=5), "cannot be greater")
  expect_error(exposure(diffnet, lags=NA))
  expect_error(exposure(diffnet, lags=c(1:2)))

})
