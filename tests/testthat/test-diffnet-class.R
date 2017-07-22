context("Diffnet class and methods")

# Checking attributes ----------------------------------------------------------
test_that("Checking attributes in new_diffnet", {
  # Generating data
  set.seed(3312)
  graph <- lapply(1990:2001, function(x) rgraph_ba(t = 9))
  attrs <- matrix(runif(10*3), ncol=3)
  dynat <- lapply(1990:2001, function(x) attrs)
  toa   <- sample(c(NA,1990:2001), 10, TRUE)

  # Dynamic attributes --------------------------------------------------------
  # Number of rows
  expect_error(
    new_diffnet(graph, toa, t0=1990, t1=2001, vertex.dyn.attrs = attrs[-1,]),
    "incorrect number of rows")

  # Length of list
  expect_error(
    new_diffnet(graph, toa, t0=1990, t1=2001, vertex.dyn.attrs = dynat[-1]),
    "equal.+slices")

  # An element has different length
  expect_error({
    dynat[[2]] <- dynat[[1]][-1,]
    dynat[[4]] <- dynat[[1]][-1,]
    new_diffnet(graph, toa, t0=1990, t1=2001, vertex.dyn.attrs = dynat)
  }, "don't have n.+rows")
  dynat <- lapply(1990:2001, function(x) attrs)

  # Elements have different classes
  expect_error({
    dynat[[6]] <- as.data.frame(dynat[[6]])
    new_diffnet(graph, toa, t0=1990, t1=2001, vertex.dyn.attrs = dynat)
  }, "should be of the same class")
  dynat[[6]] <- dynat[[1]]

  # Have different colnames
  cnames <- netdiffuseR:::make_col_names(3, TRUE)
  dynat <- lapply(dynat, function(x) {
    colnames(x) <- cnames
    x
  })
  colnames(dynat[[4]]) <- paste0(colnames(dynat[[1]]),"col")
  expect_error({
    new_diffnet(graph, toa, t0=1990, t1=2001, vertex.dyn.attrs = dynat)
  }, "have the same colname")


  # Dynamic list of vectors
  dynat <- lapply(1990:2001, function(x) attrs[,1])
  new_diffnet(graph, toa, t0=1990, t1=2001, vertex.dyn.attrs = dynat)

  # Static attributes ----------------------------------------------------------
  expect_error(
    new_diffnet(graph, toa, t0=1990, t1=2001, vertex.static.attrs = attrs[-1,]),
    "incorrect number of rows")
})

# Summary ----------------------------------------------------------------------
test_that("Summary and subsetting slices", {
  set.seed(13131)
  diffnet <- rdiffnet(100, 20)

  sum1 <- summary(diffnet, no.print = TRUE)
  sum2 <- summary(diffnet[,,7:13], no.print = TRUE)

  # Number of adopters should hold (adoption rate)
  expect_equal(sum1[nrow(sum1),c("cum_adopt")],sum2[nrow(sum2),c("cum_adopt")])
})

# diffnet.attrs errors ---------------------------------------------------------
test_that("Error messages", {
  set.seed(009)
  diffnet <- rdiffnet(80, 20, seed.nodes = "random", seed.p.adopt = .1)

  # Invalid attr.class
  expect_error(diffnet.attrs(diffnet, "vertex", attr.class = "opa"), "should only have")
  expect_error(diffnet.attrs(diffnet, "svertex", attr.class = "dyn"), "should only have")

  # Already a diffnet
  expect_message(new_diffnet(diffnet), "already.+diffnet")

})

test_that("Passing id.and.per.vars gives the right sorting", {
  # Right sort
  data(fakesurveyDyn)

  # Creating a diffnet object
  fs_diffnet <- survey_to_diffnet(
    fakesurveyDyn, "id", c("net1", "net2", "net3"), "toa", "group",
    timevar = "time", keep.isolates=TRUE, warn.coercion=FALSE)

  # Now, we extract the graph data and create a diffnet object from scratch
  graph <- fs_diffnet$graph
  ids <- fs_diffnet$meta$ids
  graph <- Map(function(g) {
    dimnames(g) <- list(ids,ids)
    g
  }, g=graph)
  attrs <- diffnet.attrs(fs_diffnet, as.df=TRUE)
  toa   <- diffnet.toa(fs_diffnet)

  # Lets apply a different sorting to the data to see if it works
  n <- nrow(attrs)
  attrs <- attrs[order(runif(n)),]

  # Now, recreating the old diffnet object (notice -id.and.per.vars- arg)
  fs_diffnet_new <- new_diffnet(graph, toa=toa, vertex.dyn.attrs=attrs,
                               id.and.per.vars = c("id", "per"))

  # Now, retrieving attributes. The 'new one' will have more (repeated)
  attrs_new <- diffnet.attrs(fs_diffnet_new, as.df=TRUE)
  attrs_old <- diffnet.attrs(fs_diffnet, as.df=TRUE)

  # Comparing elements!
  tocompare <- intersect(colnames(attrs_new), colnames(attrs_old))
  expect_true(all(attrs_new[,tocompare] == attrs_old[,tocompare], na.rm = TRUE)) # TRUE!

})

# test_that("Setting attributes", {
#   set.seed(909)
#   diffnet <- rdiffnet(80, 20, seed.nodes = "random", seed.p.adopt = .1)
#
#   expect_error(diffnet.attrs(diffnet) <- as.matrix(sample(1:4, 20, TRUE)), "different lengths")
# })

test_that("Changing toa", {
  set.seed(182321)
  diffnet <- rdiffnet(100, 10)

  # All to the first time period
  diffnet.toa(diffnet) <- 1L
  expect_output(print(diffnet), "Final prevalence\\s+[:] 1\\.0")

  # No adopters... what!?
  diffnet.toa(diffnet) <- NA
  expect_output(print(diffnet), "Final prevalence\\s+[:] 0\\.0")
})

# Checking different input classes ---------------------------------------------
test_that("new_diffnet with different graph classes", {
  # Random diffnet
  set.seed(881)
  graph0 <- rdiffnet(100, 10)

  # Getting the bits
  g_array <- as.array(graph0)
  g_dgcMa <- graph0$graph
  toa     <- diffnet.toa(graph0)
  s_attrb <- data.frame(x=graph0[["real_threshold"]])

  # Using different inputs
  dn_arr <- new_diffnet(g_array, toa=toa, vertex.static.attrs = s_attrb, t0=1, t1=10)
  dn_dgc <- new_diffnet(g_dgcMa, toa=toa, vertex.static.attrs = s_attrb, t0=1, t1=10)

  # Comparing adjmats
  test <- sapply(seq_len(10), function(x) identical(as.matrix(dn_arr$graph[[x]]), as.matrix(dn_dgc$graph[[x]])))
  expect_true(all(test))

  # Comparing the rest
  dn_arr$graph      <- dn_dgc$graph      <- NULL
  dn_arr$meta$class <- dn_dgc$meta$class <- NULL
  expect_equal(dn_arr, dn_dgc)

  # Defunct
  expect_error(diffnet.attrs(dn_arr, "hola") <- NULL, "defunct")

})

# ------------------------------------------------------------------------------
test_that("Warnings and errors", {
  set.seed(11222344)
  g <- rdiffnet(100,5)
  g[["dynamic"]] <- lapply(1:5, function(x) runif(100))

  # Messages
  expect_warning(with(g, new_diffnet(graph, as.numeric(toa))), "into integer")

  attrs <- as.matrix(g$vertex.static.attrs)
  dimnames(attrs) <- NULL

  dynattrs <- lapply(g$vertex.dyn.attrs, function(x) {
    ans <- as.matrix(x)
    dimnames(ans) <- NULL
    ans
    })
  ans <- with(g, new_diffnet(graph, toa, vertex.static.attrs = attrs,
                            vertex.dyn.attrs = dynattrs))

  # Making variables
  expect_output(print(ans), "V1 \\(1\\)")
  expect_output(print(ans), "v\\.dyn\\.1 \\(1\\)")
  expect_s3_class(ans$vertex.static.attrs, "data.frame")

  expect_error(
    with(g, new_diffnet(graph, toa, vertex.static.attrs = 1:99))
  )
  expect_error(
    with(g, new_diffnet(graph, toa, graph.attrs = vector("list",5)))
  )

  # Different method of dynamic attributes
  ans0 <- with(g, new_diffnet(graph, toa, vertex.static.attrs = attrs,
                            vertex.dyn.attrs = do.call(rbind, dynattrs)))

  ans1 <- with(g, new_diffnet(graph, toa, vertex.static.attrs = attrs,
                             vertex.dyn.attrs = do.call(rbind, dynattrs)[,1]))
  expect_equal(ans,ans0)
  expect_equal(ans,ans1)
})

# ------------------------------------------------------------------------------
test_that("graph attributes", {
  set.seed(131)
  g <- lapply(1:4, function(x) rgraph_ba(t=9))
  toa <- sample(c(NA, 1:4), 10, TRUE)

  # Less than expected
  expect_error(new_diffnet(g, toa, graph.attrs = runif(3)), "3 .+ 4")
  expect_error(new_diffnet(g, toa, graph.attrs = data.frame(runif(3))), "3 .+ 4")
  expect_output(print(
    new_diffnet(g, toa, graph.attrs = runif(4))$graph.attrs
  ), "1 .+\\n2 .+\\n3 .+\\n4")

})

