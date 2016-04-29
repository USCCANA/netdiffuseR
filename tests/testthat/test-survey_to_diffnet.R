context("Survey to diffnet")

test_that("Filling edgelist and dataset", {
  # Creating network data (1:5)
  set.seed(1231)
  n  <- 10
  el <- adjmat_to_edgelist(ring_lattice(n, 2, undirected = TRUE))

  # Creating a random data
  dat <- data.frame(id = 1:n, x = runif(n),
                    toa = as.integer(sample(c(NA, 1:5), n, TRUE)))

  # Creating diffnets
  dn <- edgelist_to_diffnet(el[,1:2], dat=dat, idvar = "id", toavar="toa")

  # Removing 2 from el
  toRemove <- which(el[,1] == 2 | el[,2] == 2)
  diffnet_filling_el <- suppressWarnings(
    edgelist_to_diffnet(el[-toRemove,1:2], dat=dat, idvar = "id",
                        toavar="toa", fill.missing = "edgelist"))

  # Removing 2 from dat
  diffnet_filling_dat <- suppressWarnings(
    edgelist_to_diffnet(el[,1:2], dat=dat[-2,], idvar = "id",
                        toavar="toa", fill.missing = "dat"))

  # Should be equal but from the graph
  diffnet0 <- dn; diffnet0$graph <- NULL
  diffnet_filling_el$graph <- NULL

  expect_equal(diffnet0, diffnet_filling_el)

  # Should be equal but from the 2 individual
  diffnet_filling_dat$graph <- NULL
  expect_equal(diffnet0[-3], diffnet_filling_dat[-3])
})

