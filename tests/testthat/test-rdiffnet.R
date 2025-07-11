context("Random diffusion graphs")

# Checking input
test_that("Input", {
  # Simple matrix --------------------------------------------------------------

  # Baseline graph
  set.seed(12312)
  x <- rgraph_ba(t=5e2-1L)

  # Using a dgCMatrix
  set.seed(131)
  x_dgCMatrix <- rdiffnet(seed.graph = x, t=10,rewire.args = list(p=c(0, rep(.1,9))))

  # Using matrix
  set.seed(131)
  x_matrix <- rdiffnet(
    seed.graph=as.array(x_dgCMatrix)[,,1],
    t=10,
    rewire.args = list(p=c(0, rep(.1,9)))
    )

  # Using a function
  set.seed(12312)
  x <- function() {
    out <- rgraph_ba(t=5e2-1L)
    set.seed(131)
    out
  }
  x_fun <- rdiffnet(seed.graph = x, t=10,rewire.args = list(p=c(0, rep(.1,9))))

  # Coercing into arrays (this is easier to compare)
  x_dgCMatrix$graph <- as.array(x_dgCMatrix)
  x_matrix$graph    <- as.array(x_matrix)
  x_fun$graph       <- as.array(x_fun)

  # Comparing
  expect_equal(x_dgCMatrix,x_matrix)
  expect_equal(x_dgCMatrix,x_fun)


  # Dynamic matrix -------------------------------------------------------------
  set.seed(12312)
  x <- rgraph_ba(t=5e2-1L)
  x <- list(x, rewire_graph(x, .1))

  set.seed(131)
  x_dgCMatrix <- rdiffnet(seed.graph = x,
                          rewire.args = list(p=0))

  set.seed(131)
  x_matrix    <- rdiffnet(seed.graph=as.array(x_dgCMatrix),
                          rewire.args = list(p=0))

  set.seed(131)
  x_diffnet   <- rdiffnet(seed.graph=x_dgCMatrix,
                          rewire.args = list(p=0))

  # Coercing into arrays (this is easier to compare)
  x_dgCMatrix$graph <- as.array(x_dgCMatrix)
  x_matrix$graph    <- as.array(x_matrix)
  x_diffnet$graph   <- as.array(x_diffnet)

  # Comparing
  expect_equal(x_dgCMatrix,x_matrix)
  expect_equal(x_dgCMatrix,x_diffnet)

})

# Seed of first adopters
test_that("All should be equal!", {
  set.seed(12131)
  g    <- rgraph_ws(20, 4, p=.3)
  set0 <- c(1,5,7,10)
  thr  <- runif(20, .4,.7)

  # Generating identical networks
  net1 <- rdiffnet(seed.graph = g, seed.nodes = set0, t = 4, rewire = FALSE,
                   threshold.dist = thr)
  net2 <- rdiffnet(seed.graph = g, seed.nodes = set0, t = 4, rewire = FALSE,
                   threshold.dist = thr)

  expect_equal(net1, net2)
})


test_that("Error and warning on rdiffnet", {

  set.seed(111)

  expect_error(rdiffnet(100, 5, threshold.dist = rep(10,10)))
  expect_error(rdiffnet(100, 5, threshold.dist = rep(10,100)), "No diffusion")
  expect_warning(rdiffnet(100, 5, threshold.dist = rep(10,100), stop.no.diff = FALSE), "No diffusion")

})

test_that("Simulation study", {

  set.seed(1)
  f <- function(x) mean(x$toa, na.rm=TRUE)
  ans0 <- suppressWarnings(rdiffnet_multiple(5, f, n=50, t=4, stop.no.diff=FALSE))

  set.seed(1)
  ans1 <- suppressWarnings(sapply(1:5, function(x) f(rdiffnet(n=50, t=4, stop.no.diff=FALSE))))

  expect_equal(ans0, ans1)

})

# Testing diffnet class across several inputs (single)
test_that("rdiffnet must run across several inputs (single)", {
  expect_s3_class(rdiffnet(100, 5), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = 0.1), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = 0.1, seed.nodes = 'random'), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.nodes = c(1, 3, 5)), "diffnet")

  # summary
  net_1 <- rdiffnet(100, 5, seed.nodes = c(1,3,5))
  expect_s3_class(summary(net_1), "data.frame")
})

# Testing diffnet class across several inputs (multiple)
test_that("rdiffnet must run across several inputs (multiple)", {
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08)), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), behavior = c('tabacco', 'alcohol')), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), seed.nodes = 'random'), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), seed.nodes = c('random', 'central')), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), threshold.dist = 0.3), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), threshold.dist = list(0.1, 0.2)), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), threshold.dist = rexp(100)), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), threshold.dist = list(rexp(100), runif(100))), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), threshold.dist = function(x) 0.3), "diffnet")
  expect_s3_class(rdiffnet(100, 5, seed.p.adopt = list(0.1, 0.08), threshold.dist = list(function(x) 0.3, function(x) 0.2)), "diffnet")

  net_2 <- rdiffnet(100, 5, seed.p.adopt = list(0.05,0.05), seed.nodes = c(1,3,5))
  expect_s3_class(summary(net_2), "data.frame")
})

test_that("All should be equal! (multiple)", {
  set.seed(12131)
  n            <- 50
  t            <- 5
  graph        <- rgraph_ws(n, 4, p=.3)
  seed.p.adopt <- list(0.1, 0.1)
  seed.nodes   <- c(1,5,7,10)
  thr          <- runif(n, .2,.4)
  thr_list     <- list(thr,thr)

  # Generating identical networks
  net1 <- rdiffnet(seed.graph = graph, seed.nodes = seed.nodes, seed.p.adopt = seed.p.adopt,
                   t = t, rewire = FALSE, threshold.dist = thr_list)

  net2 <- rdiffnet(seed.graph = graph, seed.nodes = seed.nodes, seed.p.adopt = seed.p.adopt,
                   t = t, rewire = FALSE, threshold.dist = thr_list)

  expect_equal(net1, net2)
})


test_that("toa, adopt, and cumadopt should be equal! (split_behaviors tests)", {
  set.seed(12131)
  n            <- 50
  t            <- 5
  graph        <- rgraph_ws(n, 4, p=.3)
  seed.nodes   <- c(1,5,7,10)
  thr          <- runif(n, .2,.4)

  # Generating identical networks
  net_single <- rdiffnet(seed.graph = graph, seed.nodes = seed.nodes, seed.p.adopt = 0.1,
                         t = t, rewire = FALSE, threshold.dist = thr)

  net_multiple <- rdiffnet(seed.graph = graph, seed.nodes = seed.nodes, seed.p.adopt = list(0.1, 0.1),
                           t = t, rewire = FALSE, threshold.dist = thr)

  net_single_from_multiple <- split_behaviors(net_multiple)
  net_single_from_multiple_1 <- net_single_from_multiple[[1]]

  expect_equal(net_single_from_multiple_1$toa, net_single$toa)
  expect_equal(net_single_from_multiple_1$adopt, net_single$adopt)
  expect_equal(net_single_from_multiple_1$cumadopt, net_single$cumadopt)
})

test_that("Disadoption works", {


  set.seed(1231)
  n <- 500

  d_adopt <- function(expo, cumadopt, time) {

    # Id double adopters
    ids <- which(apply(cumadopt[, time, , drop=FALSE], 1, sum) > 1)

    if (length(ids) == 0)
      return(list(integer(), integer()))

    # Otherwise, make them pick one (literally, you can only adopt
    # A single behavior, in this case, we prefer the second)
    return(list(ids, integer()))

  }

  ans_d_adopt <- rdiffnet(
    n = n, t = 10, disadopt = d_adopt,
    seed.p.adopt = list(0.1, 0.1)
    )

  tmat <- toa_mat(ans_d_adopt)
  should_be_ones_or_zeros <- tmat[[1]]$cumadopt[, 10] + tmat[[2]]$cumadopt[, 10]
  expect_true(all(should_be_ones_or_zeros %in% c(0,1)))

})
