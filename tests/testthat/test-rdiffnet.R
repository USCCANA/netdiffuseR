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
  x_matrix <- rdiffnet(seed.graph=as.array(x_dgCMatrix)[,,1], t=10, rewire.args = list(p=c(0, rep(.1,9))))

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
