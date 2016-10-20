context("Rewiring algorithms")

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

  ntimes <- 5

  # BA model
  test <- NULL
  for (i in 1:ntimes) {
    for (j in 1:ntimes) {
      graph <- rgraph_ba(t=9)
      test <- c(test, sum(graph) == sum(rewire_graph(graph, p=.5, undirected = FALSE)))
    }
  }
  expect_true(all(test))

  # Bernoulli
  test <- NULL
  for (i in 1:ntimes) {
    for (j in 1:ntimes) {
      graph  <- rgraph_er(undirected = TRUE)
      suppressWarnings(graphr <- rewire_graph(graph, p=.5, undirected = TRUE))
      test <- c(test, sum(graph) == sum(graphr))
    }
  }

  expect_true(all(test))
})

test_that("When p=1 in rewiring, Pr(j'=i) = Pr(j'=k) for all (i,k) in V", {
  # Generating seed graph
  set.seed(2991)
  n <- 1e2
  x <- ring_lattice(n, 2)

  # Simulating
  N <- 1e3
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
  n <- 5e2
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


  # # Alternating exagons (hold deg seq)
  # out <- vector(length=n)
  # for (i in 1:n) {
  #   g0 <- rgraph_ba(t = 99, self=FALSE)
  #   d0 <- dfun(g0)
  #   g1 <- netdiffuseR:::rewire_swap(g, althexagons = TRUE)
  #   d1 <- dfun(g1)
  #   out[i] <- identical(d0,d1)
  # }
  # expect_true(all(out))



})

# ------------------------------------------------------------------------------
test_that("rewire_permute", {
  set.seed(12313123)
  N <- 10
  g <- rgraph_ba(m=4, t=9)

  # Shouldn't change density
  ans <- vector("logical", N)
  for (i in 1:10)
    ans[i] <- nlinks(permute_graph(g)) == nlinks(g)
  expect_true(all(ans))

  # Shouldn't change value
  ans <- vector("logical", N)
  for (i in 1:10)
    ans[i] <- sum(permute_graph(g)) == sum(g)
  expect_true(all(ans))

  # Should be equivalent
  set.seed(1); ans0 <- permute_graph(g)
  set.seed(1); ans1 <- permute_graph(as.matrix(g))
  set.seed(1); ans2 <- permute_graph(list(g))
  set.seed(1); ans3 <- permute_graph(as.array(as.matrix(g), dim=c(10,10,1)))

  expect_equal(ans0, ans1)
  expect_equal(ans0, ans2[[1]])
  expect_equal(ans0, ans3)

  # Checking diffnet
  g <- lapply(1:5, function(x) g)
  dn <- as_diffnet(g, toa=rep(1:5, 2))

  set.seed(1); ans0 <- permute_graph(g)
  set.seed(1); ans1 <- permute_graph(dn)
  expect_equivalent(ans0, dn$graph)
})

# ------------------------------------------------------------------------------
test_that("rewire_qap", {
  set.seed(12313123)
  N <- 10
  g <- rgraph_ba(m=4, t=9)

  # Shouldn't change density
  ans <- vector("logical", N)
  for (i in 1:10)
    ans[i] <- nlinks(rewire_qap(g)) == nlinks(g)
  expect_true(all(ans))

  # Shouldn't change value
  ans <- vector("logical", N)
  for (i in 1:10)
    ans[i] <- sum(rewire_qap(g)) == sum(g)
  expect_true(all(ans))

  # Should be equivalent
  set.seed(1); ans0 <- rewire_qap(g)
  set.seed(1); ans1 <- rewire_qap(as.matrix(g))
  set.seed(1); ans2 <- rewire_qap(list(g))
  set.seed(1); ans3 <- rewire_qap(as.array(as.matrix(g), dim=c(10,10,1)))

  # Checking diffnet
  g <- lapply(1:5, function(x) g)
  dn <- as_diffnet(g, toa=rep(1:5, 2))
  dn[["dynatt"]]    <- lapply(1:5, function(x) runif(10))
  dn[["staticatt"]] <- rnorm(10)

  set.seed(1); ans0 <- rewire_qap(dn)
  set.seed(1); ans1 <- rewire_qap(dn$graph)

  expect_equivalent(ans0$graph, ans1)

  # Checking attributes ordering
  ans0 <- dn$vertex.dyn.attrs
  ans1 <- rewire_qap(dn)
  ids  <- match(nodes(dn), nodes(ans1))
  ans1 <- Map(function(x) x[ids,,drop=FALSE], ans1$vertex.dyn.attrs)
  expect_equivalent(ans0,ans1)

  ans0 <- dn$vertex.static.attrs
  ans1 <- rewire_qap(dn)
  ids  <- match(nodes(dn), nodes(ans1))
  ans1 <- ans1$vertex.static.attrs[ids,,drop=FALSE]
  expect_equivalent(ans0,ans1)

  # Toamat ordering
  ans0 <- dn$cumadopt
  ans1 <- rewire_qap(dn)
  ids  <- match(nodes(dn), nodes(ans1))
  ans1 <- ans1$cumadopt[ids,,drop=FALSE]
  expect_equivalent(ans0,ans1)

  ans0 <- dn$adopt
  ans1 <- rewire_qap(dn)
  ids  <- match(nodes(dn), nodes(ans1))
  ans1 <- ans1$adopt[ids,,drop=FALSE]
  expect_equivalent(ans0,ans1)

  ans0 <- dn$toa
  ans1 <- rewire_qap(dn)
  ids  <- match(nodes(dn), nodes(ans1))
  ans1 <- ans1$toa[ids]
  expect_equivalent(ans0,ans1)

})
