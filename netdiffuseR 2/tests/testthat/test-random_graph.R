context("Random graphs")

# ------------------------------------------------------------------------------
test_that("Bernoulli model", {
  set.seed(123)
  graph <- rgraph_er(n=3e2, p=.3)
  expect_equal(sum(graph)/length(graph), .3, tolerance=5e-3, scale=1)

  set.seed(1); ans0 <- rgraph_er(100, p=.3)
  set.seed(1); ans1 <- rgraph_er(100, p=.3, as.edgelist = TRUE)
  ans1<-edgelist_to_adjmat(ans1[,1:2])

  expect_equivalent(ans0,ans1)

  # Dynamic
  ans0 <- vector("list", 2)
  set.seed(1)
  ans0[[1]] <- rgraph_er(10)
  ans0[[2]] <- rgraph_er(10)

  set.seed(1)
  ans1 <- rgraph_er(10, 2)

  expect_equivalent(ans0,ans1)


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
  diffnet <- new_diffnet(graph, toa, t0=1, t1=3)

  # Array
  graphar <- lapply(graph, as.matrix)
  graphar <- array(unlist(graphar), dim=c(10,10,3))

  set.seed(123); x_dn <- rgraph_ba(t=2, graph=diffnet)$graph
  set.seed(123); x_ar <- rgraph_ba(t=2, graph=graphar)
  set.seed(123); x_ls <- rgraph_ba(t=2, graph=diffnet$graph)

  expect_equivalent(x_dn, x_ar)
  expect_equivalent(x_dn, x_ls)


})

# De Almeida et al. ------------------------------------------------------------
test_that("De Almeida model: Generic test", {
  set.seed(8123)
  eta <- sample(1:2, 11, TRUE)
  graph <- rgraph_ba(m0=1, t = 10, eta=eta)
  graph_big <- rgraph_ba(t = 3, graph = graph, eta=c(eta, eta[1:3]))

  expect_is(graph, "dgCMatrix")
  expect_equivalent(graph, graph_big[1:11, 1:11], info="Groing graph")

  # Checking error
  expect_error(rgraph_ba(t=10, eta=as.character(eta)), "vector")
  eta[1] <- NA
  expect_error(rgraph_ba(t=10, eta=eta), "cases")
})

test_that("De Almeida model: methods", {
  set.seed(651)
  eta <- sample(c(1,30), 10, TRUE)

  # Creating dyngraph
  graph <- vector("list", 3)
  graph[[1]] <- rgraph_ba(m0=1, t = 9, eta=eta)
  graph[[2]] <- rewire_graph(graph[[1]], .3)
  graph[[3]] <- rewire_graph(graph[[2]], .3)

  # Diffnet
  toa <- sample(c(1:3,NA),10, TRUE)
  diffnet <- new_diffnet(graph, toa, t0=1, t1=3)

  # Array
  graphar <- lapply(graph, as.matrix)
  graphar <- array(unlist(graphar), dim=c(10,10,3))

  set.seed(123); x_dn <- rgraph_ba(t=2, graph=diffnet, eta=c(eta, 1:2))$graph
  set.seed(123); x_ar <- rgraph_ba(t=2, graph=graphar, eta=c(eta, 1:2))
  set.seed(123); x_ls <- rgraph_ba(t=2, graph=diffnet$graph, eta=c(eta, 1:2))

  expect_equivalent(x_dn, x_ar)
  expect_equivalent(x_dn, x_ls)
})

# Watts-Strogatz model ---------------------------------------------------------

test_that("Watts-Strogatz model: Rewiring shouldn't change the # of elements", {
  # Generating the data
  set.seed(123)
  test1 <- NULL
  test2 <- NULL
  test3 <- NULL
  for (i in 1:50) {
    graph0 <- rgraph_ws(10,2,.1, undirected = TRUE)
    graph1 <- rgraph_ws(10,2,.5, undirected = TRUE)
    graph2 <- rgraph_ws(10,2,.5, undirected = TRUE)
    graph3 <- rgraph_ws(10,2,.8, undirected = TRUE)

    test1 <- c(test1, sum(graph0) == sum(graph1))
    test2 <- c(test2, sum(graph1) == sum(graph2))
    test3 <- c(test3, sum(graph2) == sum(graph3))
  }

  expect_true(all(test1))
  expect_true(all(test2))
  expect_true(all(test3))
})


