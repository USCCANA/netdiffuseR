context("Moran's I")

test_that("Comparing against ape::Moran.I", {
  set.seed(123)
  w <- matrix(rnorm(1000*1000), ncol=1000)
#   w <- w/rowSums(w)
#   w[!is.finite(w)] <- 0

  for (i in 1:10) {
    x <- rnorm(1000)
    expect_equal(
      moran(x,w), ape::Moran.I(x, w)$observed, scale=1, tolerance=getOption("diffnet.tol")
    )
  }
}
)

# library(ape)
# microbenchmark::microbenchmark(
#   moran(x,w),
#   Moran.I(x,w), times=100
# )
