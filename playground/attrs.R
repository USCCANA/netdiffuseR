rm(list=ls())
library(netdiffuseR)
set.seed(123)
diffnet <- rdiffnet(20, 5, seed.graph = "small-world", rewire.args = list(p=.2))


get_egonet_attrs <- function(...) UseMethod("get_egonet_attrs")

get_egonet_attrs.default <-function(graph, attrs, i=1:ncol(graph), what=c("indegree")) {
  lapply(i, function(x) {
    index <- graph[x,]
    print(index)
    attrs[index]
    })
}

get_egonet_attrs.list <- function(graph, attrs, i, t=1:length(graph), what) {
  # Checking list length
  if (any(!(t %in% 1:length(graph)))) stop("-t- must be within 1 and T.")

  lapply(graph[t], function(g) get_egonet_attrs(g, attrs, i, what))
}

get_egonet_attrs.diffnet <- function(graph, i=1:graph$meta$n, t=1:length(graph$graph), what=c("indegree")) {
  get_egonet_attrs.list(graph$graph, attrs = graph$toa, i = i, t = t, what)
}

# plot(diffnet, vertex.cex = 1, displaylabels = TRUE)
x <- get_egonet_attrs(diffnet, i=6)
