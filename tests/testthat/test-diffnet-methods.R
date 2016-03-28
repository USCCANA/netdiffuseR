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

  diffnet <- as_diffnet(graph, toa, undirected = FALSE)

  # List
  coords <- plot_diffnet(graph, toa_mat(toa)$cumadopt)
  expect_equal(dim(coords), c(11,2), info = "applying to list")

  # Array
  coords <- plot_diffnet(graphar, toa_mat(toa)$cumadopt)
  expect_equal(dim(coords), c(11,2), info = "applying to array")

  # Diffnet
  coords <- plot_diffnet(diffnet)
  expect_equal(dim(coords), c(11,2), info = "applying to diffnet")
})

################################################################################
# plot_threshold, threshold and exposure
################################################################################

test_that("Returning threshold equal to the threshold fun (plot_threshold and )", {
  # Generating a random graph
  set.seed(123)
  n <- 6
  nper <- 5

  toa <- sample(2000:(2000+nper-1), n, TRUE)
  adopt <- toa_mat(toa)

  graph <- lapply(1:nper, function(x) rgraph_ba(m0 = 1,m=1, t=n-1))
  graphar <- array(unlist(lapply(graph, as.matrix)), dim=c(n,n,nper))
  diffnet <- as_diffnet(graph, toa)

  # Computing exposure
  expos <- exposure(graph, adopt$cumadopt)
  exposar <- exposure(graphar, adopt$cumadopt)
  exposdn <- exposure(diffnet)

  # Generating graph + number
  set.seed(123)
  th   <- plot_threshold(graph, expos, toa)
  set.seed(123)
  thar <- plot_threshold(graphar, exposar, toa)
  set.seed(123)
  thdn <- plot_threshold(diffnet)


  expect_equivalent(as.matrix(th["threshold"]), threshold(expos, toa))
  expect_equivalent(as.matrix(thar["threshold"]), threshold(expos, toa))
  expect_equivalent(as.matrix(thdn["threshold"]), threshold(expos, toa))
  expect_equivalent(th, thar)
  expect_equivalent(th, thdn)
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

  diffnet <- as_diffnet(graph, toa)

  infsus <- plot_infectsuscep(graph, toa, logscale  = FALSE)
  infsusar <- plot_infectsuscep(graphar, toa, logscale  = FALSE)
  infsusdn <- plot_infectsuscep(diffnet, logscale  = FALSE)

  infect <- infection(graph, toa)
  suscep <- susceptibility(graph, toa)

  expect_equal(infsus, infsusar)
  expect_equal(infsus$infect, infect)
  expect_equal(infsus$suscept, suscep)

  expect_equal(infsusar$infect, infect)
  expect_equal(infsusar$suscept, suscep)

  expect_equal(infsusdn$infect, infect)
  expect_equal(infsusdn$suscept, suscep)

})

# Printing and summary of diffnet! ---------------------------------------------
test_that("diffnet print and summary", {
  diffnet <- lapply(1:3, rgraph_ba, m0=1,t=9)
  toa <- sample(c(2001:2003, NA), 10, TRUE)
  diffnet_und <- as_diffnet(diffnet, toa, undirected = TRUE)
  diffnet_dir <- as_diffnet(diffnet, toa, undirected = FALSE)

  expect_output(print(diffnet_und), "type.+ undirected", ignore.case=TRUE)
  expect_output(print(diffnet_dir), "type.+ directed", ignore.case=TRUE)

  expect_output(summary(diffnet_und), "Diffusion network summary")
})

test_that("summary.diffnet with slices", {
  set.seed(1313)

  for (i in 1:10) {
    net <- tryCatch(rdiffnet(30,5, seed.graph = "small-world"), error=function(e) invisible(e))
    if (inherits(net, "error")) next

    slices <- c(1,2,5)
    out1 <- summary(net, no.print=TRUE)
    out2 <- summary(net, slices=slices, no.print=TRUE)

    expect_equal(out1[slices,], out2)
  }
})

