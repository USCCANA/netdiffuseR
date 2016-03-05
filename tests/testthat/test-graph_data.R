context("Graph objects")

test_that("Checking that -graph- argument correctly passes", {

  # From the adjmat.R file -----------------------------------------------------
  # expect_error(edgelist_to_adjmat(1), "No method for")
  # expect_error(edgelist_to_adjmat("a"), "No method for")
  expect_error(adjmat_to_edgelist(1), "No method for")
  expect_error(adjmat_to_edgelist("a"), "No method for")
})

test_that("Graph class analysis", {
  # Generic error
  expect_error(classify_graph(1), "Not an object")

  # Static graphs
  expect_error(classify_graph(matrix(nrow = 0, ncol=0)), "Empty")
  expect_error(classify_graph(matrix(nrow = 2, ncol=1)), "square")
  expect_error(classify_graph(matrix(TRUE,nrow = 2, ncol=2)), "numeric")
  expect_output( print(classify_graph(matrix(1,nrow = 1, ncol=1))) , "static")

  # Dynamic graphs (simple)
  expect_error(classify_graph(array(dim = c(1,0,3))), "square")
  expect_error(classify_graph(array(dim = c(1,0,0))), "length 2")
  expect_error(classify_graph(array(TRUE,dim = c(1,1,2))), "numeric")

  # Dynamic graphs (complex)
  graph <- sapply(3:6,rgraph_ba)
  expect_error(classify_graph(graph), "dimensions .* slices .* equal")

  graph <- sapply(1:3,function(x) rgraph_ba())
  expect_output( print(classify_graph(graph)) , "dynamic")

  graph <- array(unlist(lapply(graph, as.matrix)), dim=c(11,11,3))
  expect_output( print(classify_graph(graph)) , "dynamic")

  }
)
