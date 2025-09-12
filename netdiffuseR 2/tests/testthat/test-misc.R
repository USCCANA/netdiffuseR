context("Misc functions")

# ------------------------------------------------------------------------------
test_that("The recode function for data.frame works", {
  edgelist <- data.frame(
    ego   = c(1, 10, 8),
    alter = c(10, 2, 1), stringsAsFactors = FALSE
  )

  expected_edgelist <- data.frame(
    ego   = c(1, 2, 4),
    alter = c(2, 3, 1), stringsAsFactors = FALSE
  )

  expect_equivalent(recode(edgelist), expected_edgelist)
})


# ------------------------------------------------------------------------------
test_that("as_dgCMatrix works", {
  set.seed(1231)
  x <- rgraph_er(10)

  # From matrix object
  ans_mat <- as_dgCMatrix(as.matrix(x))

  # From a network object
  ans_net <- as_dgCMatrix(network::as.network(as.matrix(x)))

  # From igraph object
  ans_ig  <- as_dgCMatrix(igraph::graph_from_adjacency_matrix(x))

  expect_equivalent(ans_mat, ans_net)
  expect_equivalent(ans_mat, ans_ig)

  expect_is(as_dgCMatrix(medInnovationsDiffNet), "list")
  expect_s4_class(as_dgCMatrix(medInnovationsDiffNet)[[1]], "dgCMatrix")

  # From array
  myarray <- array(dim=c(10,10,2))
  myarray[,,1] <- as.matrix(x)
  myarray[,,2] <- as.matrix(x)

  expect_is(as_dgCMatrix(myarray), "list")
  expect_s4_class(as_dgCMatrix(myarray)[[1]], "dgCMatrix")
})
