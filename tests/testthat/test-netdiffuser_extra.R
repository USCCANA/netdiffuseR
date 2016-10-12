context("C++ extra functions")

# ------------------------------------------------------------------------------
test_that("sp_trimatl", {
  x <- ring_lattice(10,4, undirected = TRUE)
  ans <- netdiffuseR:::sp_trimatl(x)
  expect_equal(ans + methods::getMethod("t", "dgCMatrix")(ans), x)
})

# ------------------------------------------------------------------------------
test_that("sp_diag", {
  x <- methods::as(diag(1,10,10), "dgCMatrix")
  ans1 <- netdiffuseR:::sp_diag(x,1:10)
  ans2 <- x
  Matrix::diag(ans2) <- 1:10
  expect_equal(ans1, ans2)
})

# ------------------------------------------------------------------------------
test_that("sp_as_undirected", {
  set.seed(123131312)
  ans0 <- rgraph_ba(t=9,m=3)
  ans1 <- netdiffuseR:::sp_as_undirected(ans0)
  ans0 <- ans0 + methods::getMethod("t", "dgCMatrix")(ans0)
  Matrix::diag(ans0) <- Matrix::diag(ans0)/2

  expect_equivalent(ans1, ans0)

})
