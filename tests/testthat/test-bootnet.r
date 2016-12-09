context("Network bootstrap")

# rn(list=ls())
# library(microbenchmark)
# library(netdiffuseR)
#
# set.seed(1231)
# g <- rgraph_ba(t=1e2, m=2)
# index <- sample(1:nnodes(g), nnodes(g), TRUE)
# E <- g@x
# g <- g[index,index]
#
# microbenchmark(
#   R   = netdiffuseR:::bootnet_fillselfR(g, index, E),
#   Cpp = netdiffuseR:::bootnet_fillself(g, index, E), times=1e3,
#   unit= "relative"
# )

test_that("Filling zeros", {

  set.seed(123)
  n <- 5
  g <- rgraph_ba(t = n-1, self=FALSE, m = 1)

  set.seed(1)
  ans0 <- bootnet(g, function(w, i, ...) length(w@x), R=100)

  set.seed(1)
  ans1 <- bootnet(g, function(w, i, ...) length(w@x), R=100,
                  resample.args = list(self=FALSE, useR=FALSE))

  expect_equal(ans1[-length(ans0)], ans1[-length(ans1)])
})


test_that("Methods", {
  # Generating the data
  set.seed(1291)

  # Static graphs
  graphdg <- rgraph_ba(t=9)
  graphmt <- as.matrix(graphdg)

  set.seed(123); ans0 <- resample_graph(graphdg)
  set.seed(123); ans1 <- resample_graph(graphmt)

  expect_equal(ans0, ans1)

  # Dynamic graphs
  graphls <- lapply(1:3, function(x) rgraph_ba(t=9))
  names(graphls) <- 2001:2003
  toa <- sample(c(2001:2003, NA), 10, TRUE)

  graphdn <- as_diffnet(graphls, toa, t0=2001, t1=2003)$graph
  graphar <- lapply(graphls, as.matrix)
  graphar <- array(unlist(graphar), dim=c(10,10,3),
                   dimnames = list(1:10, 1:10, 2001:2003))

  set.seed(123); ans0 <- resample_graph(graphls)
  set.seed(123); ans1 <- resample_graph(graphdn)
  set.seed(123); ans2 <- resample_graph(graphar)

  expect_equivalent(ans0, ans1)
  expect_equivalent(ans0, ans2)
})
# rn(list=ls())
# library(microbenchmark)
# library(netdiffuseR)
#
# set.seed(1231)
# g <- rgraph_ba(t=1e2, m=2)
# index <- sample(1:nnodes(g), nnodes(g), TRUE)
# E <- g@x
# g <- g[index,index]
#
# microbenchmark(
#   R   = netdiffuseR:::bootnet_fillselfR(g, index, E),
#   Cpp = netdiffuseR:::bootnet_fillself(g, index, E), times=1e3,
#   unit= "relative"
# )

rm(list=ls())
library(netdiffuseR)
set.seed(123)
n <- 5
g <- rgraph_ba(t = n-1, self=FALSE, m = 1)

set.seed(1); ans0 <- bootnet(g, function(w, ...) ifelse(inherits(w, "list"), length(w$graph@x), length(w@x)), R=100)
set.seed(1); ans1 <- bootnet(g, function(w, ...) ifelse(inherits(w, "list"), length(w$graph@x), length(w@x)), R=100,
                             resample.args = list(self=FALSE, useR=FALSE))
