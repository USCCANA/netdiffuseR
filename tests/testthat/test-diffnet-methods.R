# library(netdiffuseR)
# library(testthat)
context("plot_diffnet, threshold, and exposure")

# plot_diffnet -----------------------------------------------------------------

test_that("Should return coords of dim n x 2 (plot_diffnet)", {
  # Creating the graph
  set.seed(123)
  graph <- lapply(1:3, function(x) rgraph_ba(m0 = 1,m=1, t=10))
  toa <- sample(1:3, 11, TRUE)

  graphar <- unlist(lapply(graph, as.matrix))
  graphar <- array(graphar, dim = c(11,11,3))

  diffnet <- new_diffnet(graph, toa, undirected = FALSE)

  # List
  coords <- plot_diffnet(graph, toa_mat(toa)$cumadopt)
  expect_equal(dim(coords), c(11,2), info = "applying to list")

  # Array
  dimnames(graphar) <- list(1:nnodes(graphar), 1:nnodes(graphar))
  coords <- plot_diffnet(graphar, toa_mat(toa)$cumadopt)
  expect_equal(dim(coords), c(11,2), info = "applying to array")

  # Diffnet
  coords <- plot_diffnet(diffnet)
  expect_equal(dim(coords), c(11,2), info = "applying to diffnet")
})

test_that("More plot methods", {

  # Empty graph
  set.seed(1231)
  g <- rdiffnet(20, 4)

  ans1  <- plot(g)
  ans2 <- plot(g, layout=ans1)

  expect_equal(ans1, ans2)

  ans1 <- plot_adopters(g)
  expect_output(print(ans1), "0[.]85")

  # Invallid cex
  expect_error(plot(g, vertex.size="1"), "Invalid.+size")
  expect_error(plot_diffnet(g, vertex.size="1"), "Invalid.+size")

  expect_silent(plot_adopters(g$cumadopt))
})

# plot_threshold, threshold and exposure ---------------------------------------
context("Threshold functions")
test_that("Returning threshold equal to the threshold fun (plot_threshold and )", {
  # Generating a random graph
  set.seed(123)
  n <- 6
  nper <- 5

  toa <- sample(2000:(2000+nper-1), n, TRUE)
  adopt <- toa_mat(toa)

  graph <- lapply(1:nper, function(x) rgraph_ba(m0 = 1,m=1, t=n-1))
  graphar <- array(unlist(lapply(graph, as.matrix)), dim=c(n,n,nper))
  diffnet <- new_diffnet(graph, toa)

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

  expect_error(plot_threshold(graph), "expo.+should be pro")
  expect_error(plot_threshold(diffnet, vertex.size = "a"), "Invalid.+size")

  # Repeating cex
  expect_silent(plot_threshold(diffnet, vertex.size=.5))
  expect_warning(plot_threshold(diffnet, vertex.sides = 1.2), "integer")
  expect_error(plot_threshold(diffnet, vertex.sides = "1.2"), "integer")
  expect_error(plot_threshold(diffnet, vertex.rot = "a"), "numeric")
  expect_error(plot_threshold(diffnet, vertex.rot = rep(1,2)), "same length")
  expect_error(plot_threshold(diffnet, vertex.sides=c(1L,2L)), "same length")

})

context("Infectiousness and susceptibility (plot methods)")
# plot_infectsuscept, infection, susceptibility --------------------------------
test_that("Returning threshold equal to the infect/suscept funs", {
  # Generating a random graph
  set.seed(123)
  n <- 6
  nper <- 5
  graph <- lapply(1:nper, function(x) rgraph_ba(m0 = 1,m=1, t=n-1))
  graphar <- array(unlist(lapply(graph, as.matrix)), dim=c(n,n,nper))

  toa <- sample(2000:(2000+nper-1), n, TRUE)

  diffnet <- new_diffnet(graph, toa)

  infsus <- plot_infectsuscep(graph, toa, logscale  = FALSE, h=0)
  infsusar <- plot_infectsuscep(graphar, toa, logscale  = FALSE, h=0)
  infsusdn <- plot_infectsuscep(diffnet, logscale  = FALSE, h=0)

  infect <- infection(graph, toa)
  suscep <- susceptibility(graph, toa)

  expect_equal(infsus, infsusar)
  expect_equal(infsus$infect, infect)
  expect_equal(infsus$suscept, suscep)

  expect_equal(infsusar$infect, infect)
  expect_equal(infsusar$suscept, suscep)

  expect_equal(infsusdn$infect, infect)
  expect_equal(infsusdn$suscept, suscep)

  expect_error(plot_infectsuscep(graph), "toa.+provided")
  expect_error(suppressWarnings(plot_infectsuscep(diffnet), "undefined values"))
  data("medInnovationsDiffNet")
  expect_warning(plot_infectsuscep(medInnovationsDiffNet), "missing")


})
context("Other methods")
# Printing and summary of diffnet! ---------------------------------------------
test_that("diffnet print and summary", {
  diffnet <- lapply(1:3, rgraph_ba, m0=1,t=9)
  toa <- sample(c(2001:2003, NA), 10, TRUE)
  diffnet_und <- new_diffnet(diffnet, toa, undirected = TRUE)
  diffnet_dir <- new_diffnet(diffnet, toa, undirected = FALSE)

  expect_output(print(diffnet_und), "type.+ undirected", ignore.case=TRUE)
  expect_output(print(diffnet_dir), "type.+ directed", ignore.case=TRUE)

  expect_output(summary(diffnet_und), "Diffusion network summary")

  expect_equal(capture_output(str(diffnet)), capture_output(str(unclass(diffnet))))
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

# Concatenating diffnet --------------------------------------------------------
test_that('concatenating diffnet', {

  # Spliting and putting together
  index <- medInnovationsDiffNet[['city']]<2
  mi1 <- medInnovationsDiffNet[which(index)]
  mi2 <- medInnovationsDiffNet[which(!index)]

  mi <- c(mi1, mi2)
  test <- all(mapply(identical, mi$graph, medInnovationsDiffNet$graph))

  expect_true(test)

  # Errors
  expect_error(c(mi1, 1), 'Some objects are not of class')
  expect_error(c(mi1, mi1), 'No pair of diffnets')
  expect_error(c(mi1[,,1:4], mi2), 'same time range')
  # mi1less <- mi1[['city']] <- NULL
  # colnames(mi1less)
  # expect_error(c())

})

# Arithmetic and others---------------------------------------------------------
test_that("Arithmetic and others", {
  # Pow
  set.seed(18181)
  g <- rdiffnet(100, 3)

  ans0 <- graph_power(g, 2)
  ans1 <- g^2
  ans2 <- g
  ans2$graph <- Map(function(x) x %*% x, g$graph)

  expect_equal(lapply(ans1$graph, as.matrix), lapply(ans2$graph, as.matrix))
  expect_equal(lapply(ans0$graph, as.matrix), lapply(ans2$graph, as.matrix))

  # Substract
  ans0 <- g - c(1,2)
  ans1 <- g[-c(1,2)]
  ans2 <- g-g[-(3:100)]
  ans3 <- g - c("1","2")

  expect_equal(ans0,ans1)
  expect_equal(ans0,ans2)
  expect_equal(ans0,ans3)

  expect_error(g-"z", "right-hand side")

  # MMultiply
  ans0 <- g %*% g
  ans1 <- g
  ans1$graph <- Map(function(x) x %*% x, x=ans1$graph)

  expect_equal(ans0, ans1)

  ans0 <- g %*% g$graph[[1]]
  ans1 <- g
  ans1$graph <- Map(function(x) x %*% g$graph[[1]], x=ans1$graph)

  expect_equal(ans0, ans1)

  ans0 <- g$graph[[1]] %*% g
  ans1 <- g
  ans1$graph <- Map(function(x) g$graph[[1]] %*% x, x=ans1$graph)

  expect_equal(ans0, ans1)

  # Multiply
  ans0 <- g*2
  ans1 <- g
  ans1$graph <- Map(function(x) x*2, x=ans1$graph)

  # Multiply and transpose
  ans0 <- g*t(g)
  ans1 <- g
  ans1$graph <- Map(function(x) x*methods::getMethod("t","dgCMatrix")(x), x=ans1$graph)

  expect_equal(ans0, ans1)

  # Logical comparison
  ans0 <- g & t(g)
  ans1 <- g
  ans1$graph <- Map(function(a,b) methods::as(a & b, "dgCMatrix"), a=g$graph, b=t(g$graph))

  expect_equivalent(ans1, ans0)

  ans0 <- g | t(g)
  ans1 <- g
  ans1$graph <- Map(function(a,b) methods::as(a | b, "dgCMatrix"), a=g$graph, b=t(g$graph))

  expect_equivalent(ans1, ans0)

  # Divide by scalar
  ans0 <- g/10
  ans1 <- (1/10)/(1/g)
  expect_equivalent(ans0,ans1)

  ans0 <- diffnetLapply(g, function(cumadopt,...) mean(cumadopt))
  ans1 <- lapply(lapply(apply(g$cumadopt, 2, list), unlist), mean)

  expect_equal(ans0,ans1)
})

# nnodes and nlinks------------------------------------------------------------
test_that("nnodes and nedges", {
  set.seed(21)
  dn <- rdiffnet(20,3)

  ans0 <- nlinks(dn)
  ans1 <- nlinks(dn$graph)
  ans2 <- nlinks(as.array(dn))

  expect_equal(ans0,ans1)
  expect_equal(ans0,ans2)

  ans0 <- nnodes(dn)
  ans1 <- nnodes(dn$graph)
  ans2 <- nnodes(as.array(dn))

  expect_equal(ans0,ans1)
  expect_equal(ans0,ans2)


})
