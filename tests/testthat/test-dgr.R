context("Degree function")

# Directed graph
graph <- rbind(
  c(0,1,0),
  c(1,0,0),
  c(1,1,0)
)

oldopt <- getOption("diffnet.undirected")
options(diffnet.undirected=FALSE)

test_that("indegree of directed graph", {
  expect_equal(dgr(graph, "indegree"), cbind(c(2,2,0)))
  }
)

test_that("outdegree of directed graph", {
  expect_equal(dgr(graph, "outdegree"), cbind(c(1,1,2)))
}
)

test_that("degree of directed graph", {
  expect_equal(dgr(graph, "degree"), cbind(c(3,3,2)))
}
)

options(diffnet.undirected=TRUE)

# Undirected graph
graph <- rbind(
  c(0,1,1),
  c(1,0,1),
  c(1,1,0)
)

test_that("indegree of undirected graph", {
  expect_equal(dgr(graph, "indegree"), cbind(c(2,2,2)))
}
)

test_that("outdegree of undirected graph", {
  expect_equal(dgr(graph, "outdegree"), cbind(c(2,2,2)))
}
)

test_that("degree of undirected graph", {
  expect_equal(dgr(graph, "degree"), cbind(c(2,2,2)))
}
)

options(diffnet.undirected=oldopt)
