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
  diffnet <- as_diffnet(graph, toa, t0=1, t1=3)

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
  for (i in 1:100) {
    graph0 <- rgraph_ws(10,2,.1, undirected = TRUE)
    graph1 <- rgraph_ws(10,2,.5, undirected = TRUE)
    graph2 <- rgraph_ws(10,2,.5, undirected = TRUE)
    graph3 <- rgraph_ws(10,2,.8, undirected = TRUE)

    expect_equal(sum(graph0), sum(graph1))
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

  graphdn <- as_diffnet(graphls, toa, t0=2001, t1=2003)$graph
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
      graph  <- rgraph_er(undirected = TRUE)
      suppressWarnings(graphr <- rewire_graph(graph, p=.5, undirected = TRUE))
      expect_equal(sum(graph), sum(graphr))
    }
  }
})

test_that("When p=1 in rewiring, Pr(j'=i) = Pr(j'=k) for all (i,k) in V", {
  # Generating seed graph
  set.seed(2991)
  n <- 1e2
  x <- ring_lattice(n, 2)

  # Simulating
  N <- 1e4
  out <- lapply(seq_len(N), function(y) {
    y <- rewire_graph(x, p=1.0, self = TRUE, undirected = FALSE, both.ends = FALSE,
                      multiple = FALSE)
    y <- as.matrix(y)
    colSums(y)/sum(y)
    })

  # # Computing the probability that an j was picked.
  out <- do.call(rbind, out)
  m   <- colMeans(out)

  # Case by case (should be pretty close)
  x <- rep(0, length(m))
  names(x) <- names(m)
  # plot(m-1/n, type="l", ylim=c(-.00025,.00025))
  expect_equal(m - 1/(n), x, tolerance=.00025, check.attributes=FALSE)
})

# Rewiring degree preserve
test_that("rewire_graph_const_cpp should hold degree", {
  set.seed(18231)
  n <- 1e3
  N <- 1e2

  # Function to compute degrees
  dfun <- function(x) cbind(dgr(x, "indegree"), dgr(x, "outdegree"))

  # Directed graph
  out <- vector(length = n)
  for (i in 1:n) {
    x  <- rgraph_ba(t=N-1)
    x  <- netdiffuseR:::sp_diag(x, rep(0, N))
    d0 <- dfun(x)
    y  <- netdiffuseR:::rewire_swap(x, 100)
    d1 <- dfun(y)

    out[i] <- identical(d0, d1)
  }
  expect_equal(out, rep(TRUE, n))

  # Undirected graph
  out <- vector(length = n)
  for (i in 1:n) {
    x  <- rgraph_ws(n=N-1, k=4, p=.3)
    d0 <- dfun(x)
    y  <- netdiffuseR:::rewire_swap(x, 100, undirected = FALSE)
    d1 <- dfun(y)

    out[i] <- identical(d0, d1)
  }
  expect_equal(out, rep(TRUE, n))
})
