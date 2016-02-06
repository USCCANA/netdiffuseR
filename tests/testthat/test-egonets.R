context("Egonets")

test_that("Simple extraction throught clases", {
  set.seed(1931)
  # Generating data ----------
  n <- 100
  t <- 20
  diffnet <- rdiffnet(100,20, seed.graph = "small-world")

  # Adding attributes
  diffnet.attrs(diffnet, attr.class = "dyn") <-
    lapply(1:t, function(x) cbind(unif=runif(n), bool=runif(n) > .5))

  diffnet.attrs(diffnet, attr.class = "static") <-
    cbind(d=rowMeans(dgr(diffnet)))

  diffnet.attrs(diffnet, attr.class = "static") <-
    cbind(threshold(diffnet))

  # Computing egonets
  en_diffnet <- egonet_attrs(diffnet, fun=function(x) {
    sum(x[,"threshold"] * x[,"d"], na.rm=TRUE)/sum(x[,"d"], na.rm=TRUE)
  })

  en_list <- egonet_attrs(diffnet$graph, diffnet.attrs(diffnet), fun=function(x) {
    sum(x[,"threshold"] * x[,"d"], na.rm=TRUE)/sum(x[,"d"], na.rm=TRUE)
  })
})
