rm(list=ls())
library(foreign)

# Preparing the data -----------------------------------------------------------
# brfarmers <- read.dta("data-raw/brfarmers.dta")
load("data/brfarmers.rdata")

# Adding factors

# Subsetting
netvars <- names(brfarmers)[grepl("^net", names(brfarmers))]
othervars <- c("id", "yr", "adopt", "village")

#               storage  display     value
# variable name   type   format      label      variable label
# ----------------------------------------------------------------------------
# net31           int    %8.0g                  nomination friend 1
# net32           int    %8.0g                  nomination friend 2
# net33           int    %8.0g                  nomination friend 3
# net21           int    %8.0g                  nomination influential 1
# net22           int    %8.0g                  nomination influential 2
# net23           int    %8.0g                  nomination influential 3
# net11           int    %8.0g                  nomination practice A
# net12           int    %8.0g                  nomination practice B
# net13           int    %8.0g                  nomination practice C
# net41           int    %8.0g                  nomination coop comm proj
# yr : Year of adoption (46' -> 66')

# Influential graph
# Rogers et al. (1970)
brfarmers$yr <- as.integer(brfarmers$yr) + 1900L
brfarmers$adopt <- as.integer(brfarmers$adopt) + 1900L

# Creating an ID
brfarmers$id <- with(brfarmers, id + village*100L)
surveyed <- brfarmers$id

for (i in netvars)
  brfarmers[[i]] <- brfarmers[[i]] + brfarmers$village*100

# Removing farmes that are not part of the experiment, i.e. weren't
# surveyed.

for (i in netvars)
  brfarmers[[i]][which(!(brfarmers[[i]] %in% surveyed))] <- NA

# # Adding autoedges to farmers that are isolated
# isolated <- which(apply(brfarmers[, netvars], 1, function(x) all(is.na(x))))
# brfarmers[isolated, netvars[1]] <- brfarmers$id[isolated]

# Reshaping data
brfarmers.long <- reshape(
  brfarmers[, c(netvars, othervars)], v.names= "net",
  varying = netvars,
  timevar = "level", idvar="id", direction="long", drop = c("idold"))

library(netdiffuseR)

# Coercing the edgelist to an adjacency matrix. Here we are assuming that the
# network is constant through time.
graph <- with(
  brfarmers.long,
  edgelist_to_adjmat(cbind(id, net), t=1966-1946+1, undirected=FALSE, keep.isolates = TRUE, multiple=TRUE)
  )

# We have an extra year which we won't use
# graph <- graph[-length(graph)]

# [2016-02-24]: keep.isolates working
# # Here we are retrieving the set of individuals who actually were used in the
# # network (as these are not isolated nodes)
# used.vertex <- rownames(graph[[1]])
#
# # Create the vector (subset) of times of adoption using only the individuals
# # that are included in the adjacency matrix
# toa <- brfarmers$adopt[brfarmers$id %in% used.vertex]
toa <- brfarmers$adopt

# Creating a diffnet -----------------------------------------------------------
brfarmersDiffNet <- as_diffnet(graph, toa, vertex.static.attrs = subset(brfarmers, select=c(-id,-toa)),
                               t0=1946, t1=1966, name="Brazilian Farmers",
                               behavior="Adoption of Hybrid Corn Seeds")

# [2016-02-24]: keep.isolates working
# diffnet.attrs(brfarmersDiffNet2, "vertex", "static") <- as.matrix(subset(brfarmers, id %in% used.vertex))

# [2016-03-05]: Deprecated
# diffnet.attrs(brfarmersDiffNet, "vertex", "static") <- subset(brfarmers, select=c(-id,-toa))

# # Applying the methods
# diffnet
# summary(diffnet)
#
# d <- sqrt(dgr(diffnet$graph[[21]]))
# d <- (d - min(d) + 1)/(max(d) - min(d) + 1)*2
# plot_diffnet(diffnet, displayisolates = FALSE, displaylabels=FALSE, slices=c(1,5,9,13,17,21),
#              mai = c(0,0,0,0), vertex.cex = d)
#
# # Redoing
# graph <- with(
#   brfarmers.long,
#   edgelist_to_adjmat(cbind(id, net), undirected=FALSE, use.incomplete=FALSE, t=21)
# )
# diffnet <- as_diffnet(graph, toa, t0=1946, t1=1966)
#
# # Nice plots
# plot(diffnet, t=19)
# plot_infectsuscep(diffnet, K=5, logscale = TRUE, bins=40)
# plot_threshold(diffnet, undirected = FALSE, vertex.cex = 1/5)
# plot_adopters(diffnet)
# plot_hazard(diffnet)

save(brfarmersDiffNet, file="data/brfarmersDiffNet.rdata",
     compress = "xz")

