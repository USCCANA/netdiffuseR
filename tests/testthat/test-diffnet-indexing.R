context("Checking diffnet-indexing")

# Dynamic attributes -----------------------------------------------------------

test_that("Dynamic attributes assignment work", {
  # Preparing the data
  set.seed(13213)
  x <- rdiffnet(20, 5)

  # Calculating exposure
  expoM <- exposure(x)
  expoL <- lapply(seq_len(x$meta$nper), function(x) expoM[,x,drop=FALSE])
  expoD <- do.call(rbind, expoL)

  # Adding data
  x[["expoM"]] <- expoM
  x[["expoL"]] <- expoL
  x[["expoD"]] <- expoD

  # Checking
  eM <- x[["expoM"]]
  eL <- x[["expoL"]]
  eD <- x[["expoD"]]
  expect_equal(eM, eL)
  expect_equal(eL, eD)

  expect_error(x[["expoM"]] <- expoM[-1,], "incorrect size")
  expect_error(x[["expoM"]] <- {expoL[[2]] <- "a";expoL}, "Not all elements")

  # Removing data
  expect_warning(x[["davinci"]] <- NULL)
  ans0 <- x
  ans0[["real_threshold"]] <- NULL
  expect_error(ans0[["real_threshold"]], "No dynamic or static")

  ans0 <- x
  ans0[["expoM"]] <- NULL
  expect_error(ans0[["expoM"]], "No dynamic or static")


})

# Slices of a diffnet object ---------------------------------------------------

test_that("Getting slices", {
  set.seed(5554)
  x <- rdiffnet(50,5)

  # Adding labels
  x <- as_diffnet(x$graph, toa = (2001:2005)[diffnet.toa(x)], t0=2001, t1=2005)

  # Subslices
  s_int <- 2L:4L
  s_log <- c(FALSE, TRUE, TRUE, TRUE, FALSE)
  s_lab <- as.character(2002:2004)

  # Should be equivalent
  expect_equal(x[,,s_int], x[,,s_log])
  expect_equal(x[,,s_int], x[,,s_lab])
})

# Subsetting a diffnet object --------------------------------------------------
test_that("Subsetting diffnet objects", {
  # Generating the data
  set.seed(101)
  graph <- rdiffnet(100, 10, seed.nodes = "central")

  # Testing i subset
  i <- sample.int(100, 50)
  j <- 3:8

  subg_i <- graph[i,]
  subg_j <- graph[,,j]

  expect_equal(subg_j[i,], subg_i[,,j], info = "The order doesn't matters")
  expect_equal(graph[i,,j], subg_i[,,j], info = "The order doesn't matters 2")
  expect_error(graph[,,c(1, 6:10)], "k- must represent a range")
})


# Replacing values -------------------------------------------------------------
test_that("[<-.diffnet method", {
  set.seed(1231)
  g <- rdiffnet(20, 3)

  ans0 <- g
  ans0$graph <- lapply(ans0$graph, function(x) {x[1,1] <- 10;x})
  ans1 <- g
  ans1[1,1] <- 10

  expect_equal(ans0, ans1)

  ans0 <- g
  ans0$graph <- lapply(ans0$graph, function(x) {x[1,] <- 10;x})
  ans1 <- g
  ans1[1,] <- 10

  expect_equal(ans0, ans1)

  ans0 <- g
  ans0$graph <- lapply(ans0$graph, function(x) {x[,1] <- 10;x})
  ans1 <- g
  ans1[,1] <- 10

  expect_equal(ans0, ans1)

  ans0 <- g
  ans0$graph[[1]][,1] <- 10
  ans1 <- g
  ans1[,1,1] <- 10

  expect_equal(ans0, ans1)
})


