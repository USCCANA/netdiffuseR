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

for (i in names(EL_digraph)) {
  test_that(paste0("indegree of directed graph - ",i), {
    expect_equivalent(dgr(EL_digraph[[i]], "indegree"), cbind(c(2,2,0)))
    expect_equivalent(dgr(EL_digraph[[i]], "outdegree"), cbind(c(1,1,2)))
    expect_equivalent(dgr(EL_digraph[[i]], "degree"), cbind(c(3,3,2)))
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
    expect_equivalent(dgr(EL_digraph[[i]], "indegree"), cbind(c(2,2,2)))
    expect_equivalent(dgr(EL_digraph[[i]], "outdegree"), cbind(c(2,2,2)))
    expect_equivalent(dgr(EL_digraph[[i]], "degree"), cbind(c(2,2,2)))
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
