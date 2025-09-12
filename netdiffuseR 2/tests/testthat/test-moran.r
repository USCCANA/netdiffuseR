context("Morans I")

test_that("Get the same result as in ape", {
  set.seed(123)
  graph <- rgraph_ba(t = 199)
  w <- approx_geodesic(graph)
  x <- rnorm(200)

  # Computing Moran's I
  ans0 <- moran(x, w)

  # Comparing with the ape's package version
  ans1 <- ape::Moran.I(x, as.matrix(w))

  expect_equivalent(unclass(ans0), unclass(ans1))
})
#
# graphm <- as.matrix(graph)
# microbenchmark::microbenchmark(
#   netdiffuseR::moran(x,graph),
#   ape::Moran.I(x, graphm), times = 100,
#   unit="relative"
# )
