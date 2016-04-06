rm(list=ls())

library(sna)
library(netdiffuseR)
library(igraph)

# Loading data
data("medInnovationsDiffNet")
dat <- cbind(
  x = medInnovationsDiffNet$toa,
  y = threshold(medInnovationsDiffNet, include_censored = TRUE))

head(dat)

# Graph
graph_sna    <- medInnovationsDiffNet[,,1,drop=TRUE][[1]]
graph_igraph <- graph_from_adjacency_matrix(graph_sna)
graph_sna    <- methods::as(graph_sna, "matrix.csc")

# Some parameters
xran <- range(dat[,1], TRUE)
yran <- range(dat[,2], TRUE)
ylab <- "Threshold"
xlab <- "Time"
main <- "Time of Adoption by Network Threshold"

# Baseline graph

# SNA
gplot(graph_sna, coord = dat, xlim=xran, ylim=yran, suppress.axes = FALSE,
      xlab=xlab, ylab=ylab, main=main)

# Igraph
plot(graph_igraph, layout=dat, xlim = xran, ylim=yran, axes = TRUE,
     rescale=FALSE, edge.curved=FALSE, vertex.label=NA,
     xlab=xlab, ylab=ylab, main=main)

# netdiffuseR
plot_threshold(medInnovationsDiffNet,
               xlab=xlab, ylab=ylab, main=main, include_censored = TRUE,
               jitter.factor = c(0,0), jitter.amount = NULL)
