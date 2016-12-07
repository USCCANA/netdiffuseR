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
  ans0 <- bootnet(g, function(w, ...) length(w@x), R=100)

  set.seed(1)
  ans1 <- bootnet(g, function(w, ...) length(w@x), R=100,
                  resample.args = list(self=FALSE, useR=FALSE))

  expect_equal(ans1[-length(ans0)], ans1[-length(ans1)])
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
