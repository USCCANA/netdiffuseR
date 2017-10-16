context("Stats functions (including exposure)")

test_that("exposure calculations", {
  # Generating data
  set.seed(999)
  diffnet <- rdiffnet(40,5, seed.p.adopt = .1)

  # Default
  ans0 <- exposure(diffnet)
  ans1 <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- diffnet$graph[[x]]
     ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(1e-15+Matrix::rowSums(s))
  })))
  ans2 <- exposure(diffnet$graph, cumadopt = diffnet$cumadopt)
  ans3 <- exposure(as.array(diffnet), cumadopt = diffnet$cumadopt)

  expect_equivalent(ans0, ans1)
  expect_equivalent(ans0, ans2)
  expect_equivalent(ans0, ans3)

  # With an attribute
  X <- matrix(diffnet[["real_threshold"]], ncol=5, nrow=40, byrow = FALSE)
  ans0 <- exposure(diffnet, attrs=X)
  ans1 <- exposure(diffnet, attrs="real_threshold")
  expect_equivalent(ans0, ans1)

  expect_error(exposure(diffnet$graph, attrs="real_threshold"),"is only valid for")

  # Struct Equiv
  se <- struct_equiv(diffnet)
  se <- lapply(se, function(x) {
    ans <- methods::as(x$SE, "dgCMatrix")
    ans@x <- 1/(ans@x + 1e-20)
    ans
  })
  exp_1_diffnet <- exposure(diffnet, alt.graph = se, valued=TRUE)
  se2 <- vector("list", length(se))
  exp_1_manual <- as.matrix(do.call(cbind,lapply(diffnet$meta$pers, function(x) {
    s <- methods::as(struct_equiv(diffnet$graph[[x]])$SE, "dgCMatrix")
    s@x <- 1/(s@x + 1e-20)
    se2[[x]] <<- s
    ( s %*% diffnet$cumadopt[,x,drop=FALSE])/(Matrix::rowSums(s) +1e-20)
  })))

  expect_equivalent(unname(exp_1_diffnet), unname(exp_1_manual))

  # Lagged exposure
  ans0 <- exposure(diffnet)
  ans1 <- exposure(diffnet, lags = 1)
  ans2 <- exposure(diffnet, lags = 2)
  ans3 <- exposure(diffnet, lags = -1)

  expect_equivalent(ans0[,-5], ans1[,-1])
  expect_equivalent(ans0[,-(4:5)], ans2[,-(1:2)])
  expect_equivalent(ans0[,-1], ans3[,-5])

  expect_error(exposure(diffnet, lags=5), "cannot be greater")
  expect_error(exposure(diffnet, lags=NA))
  expect_error(exposure(diffnet, lags=c(1:2)))

})

test_that("Times of Adoption", {
  # Creating the data
  set.seed(13131)
  toa <- sample(c(NA, 1:10), 100, TRUE)

  toa_mat2 <- function(
    times, labels=NULL,
    t0=min(times, na.rm=TRUE), t1=max(times, na.rm=TRUE)) {

    # Counting number of rows
    n <- length(times)

    # Checking names
    if (length(labels)) rn <- labels
    else {
      rn <- names(times)
      if (!length(rn)) rn <- 1:n
    }

    cn <- t0:t1

    # Computing
    m <- matrix(0, nrow=n, ncol= t1-t0 + 1, dimnames = list(rn, cn))
    m[cbind(rn, times - t0 + 1)] <- 1
    m <- list(adopt=m, cumadopt=t(apply(m, 1, cumsum)))

    # Assigning names
    dimnames(m[[2]]) <- dimnames(m[[1]])

    m
  }

  expect_equal(toa_mat(toa), toa_mat2(toa))
  # library(microbenchmark)
  # tm1 <- toa_mat(toa)
  # tm2 <- toa_mat2(toa)
  #
  # microbenchmark(
  #   toa_mat(toa),
  #   toa_mat2(toa), times = 1000
  # )
})

# ------------------------------------------------------------------------------
test_that("Threshold levels", {
  set.seed(11231)
  g <- rdiffnet(n=100, t=5, seed.nodes = "central", rgraph.args=list(m=4),
                threshold.dist = function(x) .5)

  ans1 <- threshold(g)
  ans2 <- exposure(g)
  ans2 <- sapply(1:100, function(x) {
    ans2[x, g$toa[x]]
  })

  expect_equal(as.vector(ans1), ans2)
})

# ------------------------------------------------------------------------------
test_that("vertex_covariate_distance", {
  set.seed(123131231)
  n <- 20
  X <- matrix(runif(n*2, -1,1), ncol=2)
  W <- rgraph_ws(n,4,.2)

  # Mahalanobis
  D <- vertex_covariate_dist(W,X)
  D2 <- methods::as(matrix(0, n,n), "dgCMatrix")

  D2 <- methods::as(as.matrix(dist(X)), "dgCMatrix")*W

  expect_equal(sum(D2-D), 0)

  # minkowski
  D <- vertex_covariate_dist(W,X, p=1)
  D2 <- methods::as(matrix(0, n,n), "dgCMatrix")

  D2 <- methods::as(as.matrix(dist(X, method = "minkowski", p=1)), "dgCMatrix")*W

  expect_equal(sum(D2-D), 0)
})


# ------------------------------------------------------------------------------
test_that("vertex_mahalanobis_dist", {
  set.seed(123131231)
  n <- 20
  X <- matrix(runif(n*2, -1,1), ncol=2)
  G <- rgraph_ws(n,4,.2)
  W <- var(X)

  ans1 <- vertex_mahalanobis_dist(G,X, W)
  ans2 <- methods::as(matrix(0, n,n), "dgCMatrix")

  for (i in 1:n)
    for (j in 1:n)
      ans2[i,j] <- sqrt(mahalanobis(X[i,] - X[j,], FALSE, cov=W))

  ans2 <-ans2*G

  expect_equivalent(ans1,ans2)
})

# ------------------------------------------------------------------------------
test_that("vertex_covarite_compare", {
  g <- methods::as(matrix(c(0,1,1,0,0,0,0,0,0), ncol=3), "dgCMatrix")
  x <- cbind(1, 1, 3)
  expect_equal(vertex_covariate_compare(g, x, "distance")@x, 2)
  expect_equal(vertex_covariate_compare(g, x, "quaddist")@x, 4)
  expect_equal(vertex_covariate_compare(g, x, "equal")@x, 1)
  expect_equal(vertex_covariate_compare(g, x, "greater")@x, 1)
  expect_equal(vertex_covariate_compare(g, x, "greaterequal")@x, c(1,1))
  expect_equal(vertex_covariate_compare(g, x, "smaller")@x, numeric())
  expect_equal(vertex_covariate_compare(g, x, "smallerequal")@x, 1)
})

# ------------------------------------------------------------------------------
test_that("Classify adopter", {
  # Creating graph
  g <- matrix(0, ncol=3, nrow=3)
  g[cbind(1,2:3)] <- 1
  g[cbind(2:3,1)] <- 1
  g[2,3] <- 1
  g[3,2] <- 1

  g <- lapply(1:3, function(x) methods::as(g, "dgCMatrix"))

  # Cumultive adoption matrix
  toa <- 1:3

  dn <- as_diffnet(g, toa)
  ans <- ftable(classify(dn))

  expect_equal(sum(ans), 100, tolerance = .01)
  expect_equal(colSums(ans)/100, c(0, 0,1/3,1/3,1/3), tolerance = .005)
})

# ------------------------------------------------------------------------------
test_that("Approximate geodesic", {
  # RING
  g  <- ring_lattice(20, 2)
  ig <- igraph::graph_from_adjacency_matrix(g)

  # Warning
  expect_warning(approx_geodist(g, n = 1000, warn = TRUE))

  ans0 <- as.matrix(approx_geodist(g, n = 100))
  ans1 <- igraph::distances(ig, mode = "out")
  expect_equivalent(ans0,  ans1)

  # BARABASI ALBERT
  set.seed(111222)
  g  <- rgraph_ba(t=19, m =4, self = FALSE)
  ig <- igraph::graph_from_adjacency_matrix(g)

  ans0 <- approx_geodist(g)
  ans1 <- igraph::distances(ig, mode = "out")
  arenot0 <- which(as.matrix(ans0) != 0, arr.ind = TRUE)
  expect_equal(ans0[arenot0],  ans1[arenot0])

  # Watts-Strugatz
  set.seed(111222)
  g  <- rgraph_ws(n=20, k=4, p = .5)
  ig <- igraph::graph_from_adjacency_matrix(g)

  ans0 <- approx_geodist(g)
  ans1 <- igraph::distances(ig, mode = "out")
  arenot0 <- which(as.matrix(ans0) != 0, arr.ind = TRUE)
  expect_equal(ans0[arenot0],  ans1[arenot0])

})

# ------------------------------------------------------------------------------
test_that("Matrix comparison", {

  set.seed(89)
  A <- rgraph_ba(t = 9, m = 4)
  B <- rgraph_ba(t = 9, m = 4)
  A;B

  # Comparing
  ans0 <- as.matrix(matrix_compare(A,B, function(a,b) (a+b)/2))

  ans1 <- matrix(0, ncol=10, nrow=10)
  for (i in 1:10)
    for (j in 1:10)
      ans1[i,j] <- mean(c(A[i,j], B[i,j]))

  expect_equivalent(ans0[], ans1[])

})

