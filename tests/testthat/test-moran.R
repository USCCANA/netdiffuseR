context("Moran's I")

test_that("Comparing against ape::Moran.I", {
  set.seed(123)
  w <- matrix(rnorm(20*20), ncol=20)
  w <- w/rowSums(w)
  w[!is.finite(w)] <- 0

  for (i in 1:10) {
    x <- rnorm(20)
    expect_equal(
      moran(x,w), ape::Moran.I(x, w)$observed, scale=1, tolerance=getOption("diffnet.tol")
    )
  }
})
