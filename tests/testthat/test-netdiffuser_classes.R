context("Plotting functions 1 (plot_diffnet, threshold, and exposure)")

################################################################################
# plot_diffnet
################################################################################

test_that("Should return coords of dim n x 2 (plot_diffnet)", {
  # Creating the graph
  set.seed(123)
  graph <- lapply(1:3, function(x) rgraph_ba(m0 = 1,m=1, t=10))
  toa <- sample(1:3, 11, TRUE)

  graphar <- unlist(lapply(graph, as.matrix))
  graphar <- array(graphar, dim = c(11,11,3))

  # List
  coords <- plot_diffnet(graph, toa_mat(toa)$cumadopt)
  expect_equal(dim(coords), c(11,2), info = "applying to list")

  # Array
  coords <- plot_diffnet(graphar, toa_mat(toa)$cumadopt)
  expect_equal(dim(coords), c(11,2), info = "applying to array")
})

################################################################################
# plot_threshold, threshold and exposure
################################################################################

test_that("Returning threshold equal to the threshold fun (plot_threshold and )", {
  # Generating a random graph
  set.seed(123)
  n <- 6
  nper <- 5
  graph <- lapply(1:nper, function(x) rgraph_ba(m0 = 1,m=1, t=n-1))
  graphar <- array(unlist(lapply(graph, as.matrix)), dim=c(n,n,nper))

  toa <- sample(2000:(2000+nper-1), n, TRUE)
  adopt <- toa_mat(toa)

  # Computing exposure
  expos <- exposure(graph, adopt$cumadopt, undirected = FALSE)
  exposar <- exposure(graphar, adopt$cumadopt, undirected = FALSE)

  # Generating graph + number
  set.seed(123)
  th   <- plot_threshold(graph, expos, toa)
  set.seed(123)
  thar <- plot_threshold(graphar, exposar, toa)


  expect_equal(as.matrix(th["threshold"]), threshold(expos, toa))
  expect_equal(as.matrix(thar["threshold"]), threshold(expos, toa))
  expect_equal(th, thar)
})

################################################################################
# plot_infectsuscept, infection, susceptibility
################################################################################

test_that("Returning threshold equal to the infect/suscept funs", {
  # Generating a random graph
  set.seed(123)
  n <- 6
  nper <- 5
  graph <- lapply(1:nper, function(x) rgraph_ba(m0 = 1,m=1, t=n-1))
  graphar <- array(unlist(lapply(graph, as.matrix)), dim=c(n,n,nper))

  toa <- sample(2000:(2000+nper-1), n, TRUE)

  infsus <- plot_infectsuscep(graph, toa, logscale  = FALSE)
  infsusar <- plot_infectsuscep(graphar, toa, logscale  = FALSE)

  infect <- infection(graph, toa)
  suscep <- susceptibility(graph, toa)

  expect_equal(infsus, infsusar)
  expect_equal(infsus$infect, infect)
  expect_equal(infsus$suscept, suscep)
  expect_equal(infsusar$infect, infect)
  expect_equal(infsusar$suscept, suscep)

})
