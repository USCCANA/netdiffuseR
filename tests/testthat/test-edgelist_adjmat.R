context("edgelist to adjmat and vice-versa")

# Graph with attributes
edgelist <- cbind(
  c(2:4,1),
  1:4
)
w <- rowMeans(edgelist)
tim <- sample(1:4, 4, TRUE)

test_that("Length of inputs in edgelist_adjmat should match (detecting error)", {
  expect_error(edgelist_to_adjmat(edgelist, weights = w[-1]), "Error.+same length")
  expect_error(edgelist_to_adjmat(edgelist, times = tim[-1]), "Error.+same length")
})

# Undirected graph, generating the data
edgelist <- cbind(
  c(2,1,4,3),
  c(1,2,3,4)
)

edgelist_recovered <- adjmat_to_edgelist(edgelist_to_adjmat(edgelist), undirected = FALSE)
edgelist <- edgelist[order(edgelist[,1]),]
test_that("Undirected edgelist-adjmat-edgelist (shoudl hold)", {
  expect_equal(edgelist_recovered, edgelist)
})


