context("Egonets")

test_that("Simple extraction throught clases", {
  set.seed(1931)
  # Generating data ----------
  n <- 50
  t <- 10
  diffnet <- rdiffnet(50,10, seed.graph = "small-world")

  # Adding attributes

  # Dynamic
  diffnet[["unif"]]      <- runif(n*t)
  diffnet[["bool"]]      <- runif(n*t) > .5

  # Static
  diffnet[["d"]]         <- rowMeans(dgr(diffnet))
  diffnet[["threshold"]] <- threshold(diffnet)


  # Computing egonets
  en_diffnet_dn <- egonet_attrs(diffnet, fun=function(x) {
    sum(x[,"threshold"] * x[,"d"], na.rm=TRUE)/sum(x[,"d"], na.rm=TRUE)
  })

  en_diffnet_ar <- egonet_attrs(as.array(diffnet), attrs=diffnet.attrs(diffnet), fun=function(x) {
    sum(x[,"threshold"] * x[,"d"], na.rm=TRUE)/sum(x[,"d"], na.rm=TRUE)
  })

  en_list <- egonet_attrs(diffnet$graph, diffnet.attrs(diffnet), fun=function(x) {
    sum(x[,"threshold"] * x[,"d"], na.rm=TRUE)/sum(x[,"d"], na.rm=TRUE)
  })

  expect_equal(en_diffnet_dn, en_diffnet_ar)
  expect_equal(en_diffnet_dn, en_list)

  en_dgCMatrix <- egonet_attrs(
    diffnet$graph[[1]],
    diffnet.attrs(diffnet)[[1]],
    fun=function(x) {
      sum(x[,"threshold"] * x[,"d"], na.rm=TRUE)/sum(x[,"d"], na.rm=TRUE)
      })

  en_matrix <- egonet_attrs(
    as.matrix(diffnet$graph[[1]]),
    diffnet.attrs(diffnet)[[1]],
    fun=function(x) {
      sum(x[,"threshold"] * x[,"d"], na.rm=TRUE)/sum(x[,"d"], na.rm=TRUE)
      })

  expect_equal(en_dgCMatrix, en_matrix)
})

test_that("Error messages", {
  set.seed(1313131312)
  graph <- rdiffnet(100, 20)

  expect_error(egonet_attrs(graph, attrs=1), "must be a list")
  expect_error(egonet_attrs(graph, attrs=vector("list",18)), "as many elements as")
  expect_error(egonet_attrs(graph, attrs=vector("list", 20), V = -1), "out of range")
  expect_error(egonet_attrs(graph, attrs=vector("list", 20), V = NA), "is empty")
})
