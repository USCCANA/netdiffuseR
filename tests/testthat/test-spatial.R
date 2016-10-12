context("Spatial functions (beta)")

test_that("diag expansion", {
  set.seed(1231)
  dn <- rdiffnet(n=100,t=5)

  ans1 <- diag_expand(dn)
  ans2 <- diag_expand(dn$graph)

  # Checking methods
  expect_equal(ans1,ans2)

  # Checking attributes
  expect_equal(dim(ans1), c(500L,500L))
  expect_equal(nlinks(ans1), sum(unlist(nlinks(dn))))

})
