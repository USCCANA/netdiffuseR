context("Foreign function")

test_that("igraph back and forth", {
  # Checking dynamic
  set.seed(1)
  dn0 <- rdiffnet(n=100, t=4)

  ig <- diffnet_to_igraph(dn0)
  dn1 <- igraph_to_diffnet(graph.list = ig, toavar="toa")

  attrs0 <- as.data.frame(dn0)
  attrs1 <- as.data.frame(dn1)

  dn1$vertex.static.attrs <- dn1$vertex.dyn.attrs[[1]]
  dn1$vertex.dyn.attrs <- NULL
  dn0$vertex.dyn.attrs <- NULL

  expect_equal(dn0, dn1)

  expect_equal(attrs0, attrs1)

  # CHecking static
  set.seed(2)
  dn0 <- rdiffnet(n=100, t=4, rewire = FALSE)
  ig <- diffnet_to_igraph(dn0)
  dn1 <- igraph_to_diffnet(ig[[1]], toavar="toa", t0 = 1, t1 = 4)

  attrs0 <- as.data.frame(dn0)
  attrs1 <- as.data.frame(dn1)

  expect_equal(dn0, dn1)

  expect_equal(attrs0, attrs1)
})


test_that("network back and forth", {
  # Getting the data
  set.seed(1)
  dn0 <- rdiffnet(n=100, t=4)

  net <- diffnet_to_network(dn0)
  ndy <- diffnet_to_networkDynamic(dn0, netdyn.args = list(verbose=FALSE))

  dn1 <- network_to_diffnet(graph.list = net, toavar="toa")
  dn2 <- networkDynamic_to_diffnet(ndy, toavar="toa")

  dn0$graph <- as.array(dn0)
  dn1$graph <- as.array(dn1)
  dn2$graph <- as.array(dn2)

  attrs0 <- as.data.frame(dn0)
  attrs1 <- as.data.frame(dn1)
  attrs2 <- as.data.frame(dn2)

  dn0$vertex.static.attrs <- NULL
  dn0$vertex.dyn.attrs <- NULL
  dn1$vertex.static.attrs <- NULL
  dn1$vertex.dyn.attrs <- NULL
  dn2$vertex.static.attrs <- NULL
  dn2$vertex.dyn.attrs <- NULL

  expect_equal(dn0, dn1)
  expect_equal(dn0, dn2)

  expect_equal(attrs0, attrs1[,colnames(attrs0)])
  expect_equal(attrs0, attrs2[,colnames(attrs0)])
})
