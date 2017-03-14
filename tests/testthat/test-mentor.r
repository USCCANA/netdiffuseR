context("Mentor matching")

# ------------------------------------------------------------------------------
test_that("Mentor matching", {
  set.seed(1)

  g <- igraph::make_graph(~1--2,1--3,1--4,4--5,6)
  g <- igraph::as_adj(g)

  ans <- mentor_matching(g, 2)
  expect_equal(sort(unique(ans$match)), as.character(c(1,4)))
})


# ------------------------------------------------------------------------------
test_that("Plot mentor", {
  set.seed(1)

  g <- igraph::make_graph(~1--2,1--3,1--4,4--5,6)
  g <- igraph::as_adj(g)
  ans <- mentor_matching(g, 2)

  expect_silent(plot(ans))
})
