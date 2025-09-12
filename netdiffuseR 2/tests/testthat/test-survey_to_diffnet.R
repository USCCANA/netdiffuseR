context("Survey to diffnet")

# ------------------------------------------------------------------------------
test_that("Error messages", {
  # Class checkings
  expect_error(
    netdiffuseR:::check_var_class_and_coerce("Murder", USArrests, "character", "integer"),
    "is not supported"
  )

  # Class cohercing
  expect_warning(
    netdiffuseR:::check_var_class_and_coerce("Murder", USArrests, "numeric", "integer",
                                             TRUE)
  )

  # Can't find variables
  expect_error(
    edgelist_to_diffnet(1, dat=USArrests, idvar="a", toavar="b", timevar = "c"),
    "can\'t be found"
  )

  expect_error(
    survey_to_diffnet(fakesurvey, "id1", c("net1", "net2", "net3"), "toa",
                      "group"),
    "can\'t be found"
  )

  # Infinite toavar
  dat <- fakesurvey
  dat$toa[1] <- Inf
  expect_error(
    survey_to_diffnet(dat, "id", c("net1", "net2", "net3"), "toa",
                      "group"),
    "Inf values"
  )

  dat <- fakesurvey
  dat$id <- with(dat, group*100 + id)
  dat$toa[1] <- Inf
  expect_error(
  suppressWarnings(
    edgelist_to_diffnet(fakeEdgelist[,1:2], dat=dat, toavar = "toa", idvar="id")
  ), "Inf values found"
  )

  dat <- fakesurveyDyn
  dat$time[1] <- Inf
  expect_error(
    survey_to_diffnet(dat, "id", c("net1", "net2", "net3"), "toa",
                      "group", timevar = "time"),
    "Inf values"
  )

  # Can't have NAs in timevar
  dat <- fakesurveyDyn
  dat$time[1] <- NA
  expect_error(
    survey_to_diffnet(dat, "id", c("net1", "net2", "net3"), "toa",
                      "group", timevar = "time", warn.coercion = FALSE),
    "timevar- have missing"
  )

  # Time invariant TOA
  dat <- fakesurveyDyn
  dat$toa[10] <- 1990
  expect_error(
    survey_to_diffnet(dat, "id", c("net1", "net2", "net3"), "toa",
                      "group", timevar = "time", warn.coercion = FALSE),
    "toavar- is not time-invariant"
  )

  # Times out of range
  dat <- fakesurveyDyn
  dat$id <- with(dat, group*100 + id)
  dat$toa[1] <- 2001
  expect_error(
    suppressWarnings(
      edgelist_to_diffnet(fakeEdgelist[,1:2], dat=dat, toavar = "toa", idvar="id",
                          timevar = "time")
    ), "Invalid range"
  )

  dat <- fakesurveyDyn
  dat$id <- with(dat, group*100 + id)
  expect_error(
    suppressWarnings(
      edgelist_to_diffnet(fakeEdgelist[,1:2], t0=rep(2001,nrow(fakeEdgelist)),dat=dat, toavar = "toa", idvar="id",
                          timevar = "time")
    ), "Time ranges in .*should be the same"
  )

  # There can't be NA in the ids
  dat <- fakesurveyDyn
  dat$id[c(1,5)] <- NA
  expect_error(
    edgelist_to_diffnet(fakeEdgelist[,1:2], dat=dat,
                        idvar="id", toavar="toa", timevar="time"),
    "1, 5"
  )

  dat <- fakesurveyDyn
  expect_error(
    edgelist_to_diffnet(as.matrix(fakeEdgelist[,1:2]), dat=dat, warn.coercion = FALSE,
                        idvar="id", toavar="toa", timevar="time", fill.missing = TRUE),
    "currently supported"
  )

  dat <- fakesurveyDyn
  dat$time[1] <- Inf
  expect_error(
    edgelist_to_diffnet(as.matrix(fakeEdgelist[,1:2]), dat=dat, warn.coercion = FALSE,
                        idvar="id", toavar="toa", timevar="time", fill.missing = "both"),
    "Inf"
  )
})

# ------------------------------------------------------------------------------
test_that("Reproducing observed", {

  # What are the netvars
  netvars <- names(medInnovations)[grepl("^net", names(medInnovations))]

  medInnovationsDiffNet2 <- survey_to_diffnet(
    medInnovations,
    "id", netvars, "toa", "city",
    warn.coercion=FALSE, multiple = TRUE)

  # Graph and attributes
  expect_equal(medInnovationsDiffNet$graph, medInnovationsDiffNet2$graph)
  expect_equal(diffnet.toa(medInnovationsDiffNet), diffnet.toa(medInnovationsDiffNet2))
  expect_equal(diffnet.attrs(medInnovationsDiffNet),
               diffnet.attrs(medInnovationsDiffNet2))

})

# ------------------------------------------------------------------------------
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
  diffnet0 <- dn # diffnet0$graph <- NULL
  # diffnet_filling_el$graph <- NULL
  diffnet_filling_el$graph <- diffnet0$graph

  expect_equal(as.array(diffnet0), as.array(diffnet_filling_el))

  # Should be equal but from the 2 individual
  # diffnet_filling_dat$graph <- NULL
  expect_equal(diffnet0[-3], diffnet_filling_dat[-3])
})

