context("Exposure")

test_that("exposure calculations", {
  # Generating data
  set.seed(999)
  diffnet <- rdiffnet(40,5, seed.p.adopt = .1)

  # Default
  exp_0_diffnet <- exposure(diffnet)
  exp_0_manual <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- diffnet$graph[[x]]
     ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(1e-15+Matrix::rowSums(s))
  })))

  expect_equivalent(exp_0_diffnet, exp_0_manual)

  # Struct Equiv
  se <- struct_equiv(diffnet)
  se <- lapply(se, function(x) methods::as((x$SE)^(-1), "dgCMatrix"))
  exp_1_diffnet <- exposure(diffnet, alt.graph = se)
  exp_1_manual <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- struct_equiv(diffnet$graph[[x]])$SE^(-1)
    s[!is.finite(s)] <- 0
    ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(1e-15+base::rowSums(s))
  })))

  # expect_equivalent(exp_1_diffnet, exp_1_manual)
})

# test_that("Times of Adoption", {
#   # Creating the data
#   set.seed(13131)
#   toa <- sample(c(NA, 1:10), 100, TRUE)
#
#   toa_mat2 <- function(
#     times, labels=NULL,
#     t0=min(times, na.rm=TRUE), t1=max(times, na.rm=TRUE)) {
#
#     # Counting number of rows
#     n <- length(times)
#
#     # Checking names
#     if (length(labels)) rn <- labels
#     else {
#       rn <- names(times)
#       if (!length(rn)) rn <- 1:n
#     }
#
#     cn <- t0:t1
#
#     # Computing
#     m <- matrix(0, nrow=n, ncol= t1-t0 + 1, dimnames = list(rn, cn))
#     m[cbind(cn, times - t0 + 1)] <- 1
#     m <- list(adopt=m, cumadopt=t(apply(m, 1, cumsum)))
#
#     # Assigning names
#     dimnames(m[[2]]) <- dimnames(m[[1]])
#
#     m
#   }
#
#   expect_equal(toa_mat(toa), toa_mat2(toa))
# #   library(microbenchmark)
# #   tm1 <- toa_mat(toa)
# #   tm2 <- toa_mat2(toa)
# #
# #   microbenchmark(
# #     toa_mat(toa),
# #     toa_mat2(toa), times = 1
# #   )
# })
