context("Degree function")

################################################################################
# Static graphs
################################################################################

# Directed graph ---------------------------------------------------------------
graph <- rbind(
  c(0,1,0),
  c(1,0,0),
  c(1,1,0)
)
dimnames(graph) <- list(letters[1:3],letters[1:3])

EL_digraph <- list(
  `matrix` = graph,
  `dgCMatrix` = as(graph, "dgCMatrix"))

oldopt <- getOption("diffnet.undirected")
options(diffnet.undirected=FALSE)

# Function to create diffnet_degSeq, matrix
asdegseq <- function(x) {
  structure(cbind(x), class=c("diffnet_degSeq", "matrix"))
}

for (i in names(EL_digraph)) {
  test_that(paste0("indegree of directed graph - ",i), {

    expect_equivalent(dgr(EL_digraph[[i]], "indegree"), asdegseq(c(2,2,0)))
    expect_equivalent(dgr(EL_digraph[[i]], "outdegree"), asdegseq(c(1,1,2)))
    expect_equivalent(dgr(EL_digraph[[i]], "degree"), asdegseq(c(3,3,2)))
  }
  )
}

options(diffnet.undirected=TRUE)

# Undirected graph -------------------------------------------------------------
graph <- rbind(
  c(0,1,1),
  c(1,0,1),
  c(1,1,0)
)

dimnames(graph) <- list(letters[1:3],letters[1:3])

EL_digraph <- list(`matrix` = graph, `dgCMatrix` = as(graph, "dgCMatrix"))

for (i in names(EL_digraph)) {
  test_that(paste0("Degree of undirected graph - ",i), {
    expect_equivalent(dgr(EL_digraph[[i]], "indegree"), asdegseq(c(2,2,2)))
    expect_equivalent(dgr(EL_digraph[[i]], "outdegree"), asdegseq(c(2,2,2)))
    expect_equivalent(dgr(EL_digraph[[i]], "degree"), asdegseq(c(2,2,2)))
  }
  )
}
options(diffnet.undirected=oldopt)

################################################################################
# Dynamic graphs
################################################################################

# Directed graph ---------------------------------------------------------------
graph <- rbind(
  c(0,1,0),
  c(1,0,0),
  c(1,1,0)
)
dimnames(graph) <- list(letters[1:3],letters[1:3])

listgraph <- as(graph, "dgCMatrix")
listgraph <- list(`2001`=listgraph, `2002`=listgraph, `2003`=listgraph)
arraygraph <- array(graph, dim = c(dim(graph), 3),
               dimnames=list(rownames(graph), colnames(graph), 2001:2003))

toa <- c(2001L,2001L, 2003L)

EL_digraph <- list(`array` = arraygraph, `list` = listgraph,
                   `diffnet`=as_diffnet(listgraph, toa))

oldopt <- getOption("diffnet.undirected")
options(diffnet.undirected=FALSE)

# Comparing outputs for different classes
test_that("Either as an array or as a list, dgr should work (directed)", {
  expect_identical(dgr(EL_digraph[[1]]), dgr(EL_digraph[[2]]))
})

options(diffnet.undirected=oldopt)

options(diffnet.undirected=TRUE)

# Undirected graph -------------------------------------------------------------
graph <- rbind(
  c(0,1,1),
  c(1,0,1),
  c(1,1,0)
)

dimnames(graph) <- list(letters[1:3],letters[1:3])
toa <- c(2001L,2001L, 2003L)

listgraph <- as(graph, "dgCMatrix")
listgraph <- list(`2001`=listgraph, `2002`=listgraph, `2003`=listgraph)
arraygraph <- array(graph, dim = c(dim(graph), 3),
                    dimnames=list(rownames(graph), colnames(graph), 2001:2003))

EL_digraph <- list(`array` = arraygraph, `list` = listgraph,
                   `diffnet`=as_diffnet(listgraph, toa))

# Comparing outputs for different classes
test_that("Either as an array or as a list, dgr should work (directed)", {
  expect_identical(dgr(EL_digraph[[1]]), dgr(EL_digraph[[2]]))
})

options(diffnet.undirected=oldopt)

# igraph and statnet -----------------------------------------------------------
test_that("dgr on igraph and statnet", {
  set.seed(1313)
  g <- rgraph_ba(t=9, m=2, self=FALSE)

  ans0 <- dgr(g)
  ans1 <- suppressWarnings(dgr(igraph::graph_from_adjacency_matrix(g)))
  ans2 <- dgr(network::network(as.matrix(g)))

  expect_equal(ans0,ans1)
  expect_equal(ans0,ans2)

  expect_silent(plot(ans0))

})
