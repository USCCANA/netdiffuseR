context("Random graphs")

test_that("Bernoulli model", {
  set.seed(123)
  graph <- rgraph_er(n=5e2, p=.3)
  expect_equal(sum(graph)/length(graph), .3, tolerance=5e-4, scale=1)
})

# Barabasi-albert (need more sophisticated test) -------------------------------
test_that("Barabasi-Albert model: generic test", {
  set.seed(8123)
  graph <- rgraph_ba(m0=1, t = 10)
  graph_big <- rgraph_ba(t = 3, graph = graph)

  expect_is(graph, "dgCMatrix")
  expect_equivalent(graph, graph_big[1:11, 1:11], info="Groing graph")
})

test_that("Barabasi-Albert model: methods", {
  set.seed(651)

  # Creating dyngraph
  graph <- vector("list", 3)
  graph[[1]] <- rgraph_ba(m0=1, t = 9)
  graph[[2]] <- rewire_graph(graph[[1]], .3)
  graph[[3]] <- rewire_graph(graph[[2]], .3)

  # Diffnet
  toa <- sample(c(1:3,NA),10, TRUE)
  diffnet <- as_diffnet(graph, toa)

  # Array
  graphar <- lapply(graph, as.matrix)
  graphar <- array(unlist(graphar), dim=c(10,10,3))

  set.seed(123); x_dn <- rgraph_ba(t=2, graph=diffnet)$graph
  set.seed(123); x_ar <- rgraph_ba(t=2, graph=graphar)

  expect_equivalent(x_dn, x_ar)


})

# Watts-Strogatz model ---------------------------------------------------------

test_that("Watts-Strogatz model: Rewiring shouldn't change the # of elements", {
  # Generating the data
  set.seed(123)
  for (i in 1:20) {
    graph1 <- rgraph_ws(10,2,.5)
    graph2 <- rgraph_ws(10,2,.5)
    graph3 <- rgraph_ws(10,2,.8)

    expect_equal(sum(graph1), sum(graph2))
    expect_equal(sum(graph2), sum(graph3))
  }
})


# Rewiring ---------------------------------------------------------------------
test_that("Rewiring methods", {
  # Generating the data
  set.seed(1291)

  # Static graphs
  graphdg <- rgraph_ba(t=9)
  graphmt <- as.matrix(graphdg)

  set.seed(123); graphdg <- rewire_graph(graphdg, .3)
  set.seed(123); graphmt <- rewire_graph(graphmt, .3)

  expect_equal(graphdg, graphmt)

  # Dynamic graphs
  graphls <- lapply(1:3, function(x) rgraph_ba(t=9))
  names(graphls) <- 2001:2003
  toa <- sample(c(2001:2003, NA), 10, TRUE)

  graphdn <- as_diffnet(graphls, toa)$graph
  graphar <- lapply(graphls, as.matrix)
  graphar <- array(unlist(graphar), dim=c(10,10,3),
                   dimnames = list(1:10, 1:10, 2001:2003))

  set.seed(123); graphls <- rewire_graph(graphls, .3)
  set.seed(123); graphdn <- rewire_graph(graphdn, .3)
  set.seed(123); graphar <- rewire_graph(graphar, .3)

  expect_equal(graphls, graphdn)
  expect_equal(graphar, graphdn)


})

test_that("Rewiring must hold graph's density", {
  set.seed(1231)

  ntimes <- 10

  # BA model
  for (i in 1:ntimes) {
    for (j in 1:ntimes) {
      graph <- rgraph_ba(t=9)
      expect_equal(sum(graph), sum(rewire_graph(graph, p=.5, undirected = FALSE)))
    }
  }

  # Bernoulli
  for (i in 1:ntimes) {
    for (j in 1:ntimes) {
      graph <- rgraph_er(undirected = TRUE)
      expect_equal(sum(graph), sum(rewire_graph(graph, p=.5, undirected = TRUE)))
    }
  }
})
