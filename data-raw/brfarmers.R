rm(list=ls())
library(foreign)

# Preparing the data -----------------------------------------------------------
brfarmers <- read.dta("data-raw/brfarmers.dta")

# Subsetting
netvars <- names(brfarmers)[grepl("^net", names(brfarmers))]
othervars <- c("id", "yr", "adopt", "village")
brfarmers <- brfarmers[, c(netvars, othervars)]

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

# Reshaping data
brfarmers.long <- reshape(
  brfarmers, v.names= "net",
  varying = netvars,
  timevar = "level", idvar="id", direction="long", drop = c("idold"))

library(netdiffuseR)

# Coersing the edgelist to an adjacency matrix. Here we are assuming that the
# network is constant through time.
graph <- with(
  brfarmers.long,
  edgelist_to_adjmat(cbind(id, net), time=yr, undirected=FALSE, use.incomplete=FALSE)
  )

# We have an extra year which we won't use
# graph <- graph[-length(graph)]

# Here we are retrieving the set of individuals who actually were used in the
# network (as these are not isolated nodes)
used.vertex <- rownames(graph[[1]])

# Create the vector (subset) of times of adoption using only the individuals
# that are included in the adjacency matrix
toa <- brfarmers$adopt[brfarmers$id %in% used.vertex]

# Creating a diffnet -----------------------------------------------------------
diffnet <- as_diffnet(graph, toa)

# Applying the methods
diffnet
summary(diffnet)

plot_diffnet(diffnet, displayisolates = FALSE, displaylabels=FALSE, slices=c(1,2,11,10,15,16,19,20),
             mai = c(0,0,0,0), vertex.cex = "indegree")

# Redoing
graph <- with(
  brfarmers.long,
  edgelist_to_adjmat(cbind(id, net), undirected=FALSE, use.incomplete=FALSE, t=20)
)
diffnet <- as_diffnet(graph, toa)

# Nice plots
plot_infectsuscep(diffnet, K=5, logscale = TRUE, bins=40)
plot_threshold(diffnet, undirected = FALSE, vertex.cex = 1/5)
plot_adopters(diffnet)
plot_hazard(diffnet)

