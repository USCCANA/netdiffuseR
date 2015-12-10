context("Random graphs")

test_that("Bernoulli random graphs", {
  set.seed(123)
  graph <- rgraph_er(n=5e2, p=.3)
  expect_equal(sum(graph)/length(graph), .3, tolerance=5e-4, scale=1)
})

# Barabasi-albert (need more sophisticated test) -------------------------------
test_that("BA random graphs", {
  graph <- rgraph_ba(m0=1, t = 10)
  expect_is(graph, "dgCMatrix")
})
