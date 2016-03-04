rm(list=ls())

library(netdiffuseR)
set.seed(111)

graph <- rdiffnet(10,5)

`[[.diffnet` <- function(x, name) {
    x$vertex.static.attrs[[name]]
}

`[[<-.diffnet` <- function(x, i, j, value) {
  x$vertex.static.attrs[[i]][j] <- value
  x
}

`[.diffnet` <- function(x, i, j, k) {
  # Checking ids
  if (missing(i)) i <- seq_len(x$meta$n)
  if (missing(j)) j <- seq_len(x$meta$n)
  if (missing(k)) k <- seq_len(x$meta$nper)

  lapply(x$graph[k], "[", i=i, j=j, drop=FALSE)
}

`[<-.diffnet` <- function(x, i, j, k, value) {
  # Checking ids
  if (missing(i)) i <- seq_len(x$meta$n)
  if (missing(j)) j <- seq_len(x$meta$n)
  if (missing(k)) k <- seq_len(x$meta$nper)

  for (l in k)
    x$graph[[l]][i,j] <- value

  x
}


graph[["real_threshold"]]
graph[1,,1]
graph[,,1]

graph[1,,] <- -5
