rm(list=ls())
library(foreign)

# Preparing the data -----------------------------------------------------------
medInnovations <- read.dta("data-raw/mi_v2.dta")

# Creating unique ids (including for the network data)
othervars <- c("id", "toa", "city")
netvars <- names(medInnovations)[grepl("^net", names(medInnovations))]
for (i in c("id", netvars))
  medInnovations[[i]] <- medInnovations[[i]] + medInnovations$city*1000

# Leaving unsurveyed individuals with NA
surveyed <- medInnovations$id
for (i in netvars)
  medInnovations[[i]][which(!(medInnovations[[i]] %in% surveyed))] <- NA

# [2016-02-24]: keep.isolates working
# # Adding autoedges to farmers that are isolated, we need to do this otherwize
# # these will be dropped when calling the function -edgelist_to_adjmat-. Notice
# # that this does not imply that the graph will have autoedges. (see manual)
# isolated <- which(apply(medInnovations[, netvars], 1, function(x) all(is.na(x))))
# medInnovations[isolated, netvars[1]] <- medInnovations$id[isolated]

# Reshaping data
medInnovations.long <- reshape(
  medInnovations[,c(othervars, netvars)], v.names= "net",
  varying = netvars,
  timevar = "level", idvar="id", direction="long")


library(netdiffuseR)

# Coercing the edgelist to an adjacency matrix. Here we are assuming that the
# network is constant through time.
graph <- with(
  medInnovations.long,
  edgelist_to_adjmat(cbind(id, net), t=18,undirected=FALSE, keep.isolates = TRUE, multiple=TRUE)
)

# [2016-02-24]: keep.isolates working
# # Here we are retrieving the set of individuals who actually were used in the
# # network (as these are not isolated nodes)
# used.vertex <- rownames(graph[[1]])
#
# # Create the vector (subset) of times of adoption using only the individuals
# # that are included in the adjacency matrix
# toa <- medInnovations$toa[medInnovations$id %in% used.vertex]
toa <- medInnovations$toa

# Creating a diffnet -----------------------------------------------------------
medInnovationsDiffNet <- as_diffnet(graph, toa,
                                    vertex.static.attrs = subset(medInnovations, select=c(-id,-toa)))
# [2016-02-24]: keep.isolates working
# diffnet.attrs(medInnovationsDiffNet) <-
#   subset(medInnovations, id %in% used.vertex)

# [2016-03-05]: Deprecated
# diffnet.attrs(medInnovationsDiffNet) <-
#   subset(medInnovations, select=c(-id,-toa))

# # Applying the methods
# diffnet
# summary(diffnet)
#
# d <- sqrt(dgr(diffnet$graph[[diffnet$meta$nper]]))
# d <- (d - min(d) + 1)/(max(d) - min(d) + 1)*2
# plot_diffnet(diffnet, displayisolates = FALSE, displaylabels=FALSE,
#              slices=seq(1, diffnet$meta$nper, length.out = 6),
#              mai = c(0,0,0,0), vertex.cex = d)
#
# # Nice plots
# plot(diffnet, t=18)
# plot_infectsuscep(diffnet, K=5, logscale = TRUE, bins=40)
# plot_threshold(diffnet, undirected = FALSE, vertex.cex = 1/5)
# plot_adopters(diffnet)
# plot_hazard(diffnet)

save(medInnovationsDiffNet, file="data/medInnovationsDiffNet.rdata",
     compress = "xz")

# z<-struct_test(
#   medInnovationsDiffNet,
#   function(x) mean(threshold(x), na.rm = TRUE),
#   R=1000, ncpus=8, parallel="multicore")
