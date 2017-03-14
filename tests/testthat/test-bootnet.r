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

# ------------------------------------------------------------------------------
test_that("Filling zeros", {

  set.seed(123)
  n <- 5
  g <- rgraph_ba(t = n-1, self=FALSE, m = 1)

  set.seed(1)
  ans0 <- bootnet(g, function(w, i, ...) length(w@x), R=100,
                  resample.args = list(self=FALSE, useR=FALSE))

  set.seed(1)
  ans1 <- bootnet(g, function(w, i, ...) length(w@x), R=100,
                  resample.args = list(self=FALSE, useR=TRUE))

  expect_equal(ans1[-length(ans0)], ans1[-length(ans1)])
})

# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
test_that("diffnet_bootnet methods", {
  set.seed(1222)
  x <- rgraph_ba(t=19, m=1)
  ans <- bootnet(x, function(g,...) mean(dgr(g)), R=50)

  expect_output(print(ans), "Network Bootstrap")
  expect_silent(hist(ans, ask=FALSE))
  expect_s3_class(c(ans, ans), "diffnet_bootnet")
  expect_output(print(c(ans, ans)), ": 100")

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

# rm(list=ls())
# library(netdiffuseR)
# set.seed(123)
# n <- 5
# g <- rgraph_ba(t = n-1, self=FALSE, m = 1)
#
# set.seed(1); ans0 <- bootnet(g, function(w, ...) ifelse(inherits(w, "list"), length(w$graph@x), length(w@x)), R=100)
# set.seed(1); ans1 <- bootnet(g, function(w, ...) ifelse(inherits(w, "list"), length(w$graph@x), length(w@x)), R=100,
#                              resample.args = list(self=FALSE, useR=FALSE))
#

# library(netdiffuseR)
# n <- 100
# G <- rgraph_ws(n=n, p = .5, undirected = FALSE, self=FALSE, k=4)
# G <- list(G,G)
# Y1 <- sample(c(0,1), n, TRUE)
# Y2 <- Y1
# Y2[Y2==0] <- sample(c(0,1), sum(Y2==0), TRUE)
# X <- runif(n*2)
# dat <- data.frame(Y=c(Y1, Y2), X, year=c(rep(1,n), rep(2,n)))
# dat1 <- subset(dat, year==1)
# dat2 <- subset(dat, year==2)
#
# ans1 <- bootnet(G, function(g, idx) {
#   d <- rbind(dat1[idx,], dat2[idx,])
#   suppressMessages(netmatch(d, g, "year", "Y", "X", treat_thr = 3, method="cem")$fATT)
# }, R=1e3)
#
# ans2 <- struct_test(G, function(g, idx) {
#   suppressMessages(netmatch(dat, g, "year", "Y", "X", treat_thr = 3, method="cem")$fATT)
# }, R=1e3)
#
