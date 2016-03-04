context("Egonets")

test_that("Simple extraction throught clases", {
  set.seed(1931)
  # Generating data ----------
  n <- 100
  t <- 20
  diffnet <- rdiffnet(100,20, seed.graph = "small-world")

  # Adding attributes
  diffnet.attrs(diffnet, attr.class = "dyn") <-
    lapply(1:t, function(x) cbind(unif=runif(n)))
  diffnet.attrs(diffnet, attr.class = "dyn") <-
    lapply(1:t, function(x) cbind(bool=runif(n) > .5))

  diffnet.attrs(diffnet, attr.class = "static") <-
    cbind(d=rowMeans(dgr(diffnet)))

  diffnet.attrs(diffnet, attr.class = "static") <-
    cbind(threshold(diffnet))

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
})

test_that("Error messages", {
  set.seed(1313131312)
  graph <- rdiffnet(100, 20)

  expect_error(egonet_attrs(graph, attrs=1), "must be a list")
  expect_error(egonet_attrs(graph, attrs=vector("list",18)), "as many elements as")
})
