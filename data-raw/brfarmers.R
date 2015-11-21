library(foreign)

x <- read.dta("data-raw/brfarmers.dta")

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
brfarmers <- subset(x, select=c(
  idold, net31, net32, net33, net21, net22, net23, net11, net12, net13, yr, village))
brfarmers$yr <- as.integer(brfarmers$yr) + 1900L

# Creating an ID
surveyed        <- with(brfarmers, idold + village*100L)
brfarmers$id    <- surveyed
brfarmers$net31 <- with(brfarmers, net31 + village*100L)
brfarmers$net32 <- with(brfarmers, net32 + village*100L)
brfarmers$net33 <- with(brfarmers, net33 + village*100L)
brfarmers$net21 <- with(brfarmers, net21 + village*100L)
brfarmers$net22 <- with(brfarmers, net22 + village*100L)
brfarmers$net23 <- with(brfarmers, net23 + village*100L)
brfarmers$net11 <- with(brfarmers, net11 + village*100L)
brfarmers$net12 <- with(brfarmers, net12 + village*100L)
brfarmers$net13 <- with(brfarmers, net13 + village*100L)

# Removing farmes that are not part of the experiment
brfarmers$net31[which(!(brfarmers$net31 %in% surveyed))]  <- NA
brfarmers$net32[which(!(brfarmers$net32 %in% surveyed))]  <- NA
brfarmers$net33[which(!(brfarmers$net33 %in% surveyed))]  <- NA
brfarmers$net21[which(!(brfarmers$net21 %in% surveyed))]  <- NA
brfarmers$net22[which(!(brfarmers$net22 %in% surveyed))]  <- NA
brfarmers$net23[which(!(brfarmers$net23 %in% surveyed))]  <- NA
brfarmers$net11[which(!(brfarmers$net11 %in% surveyed))]  <- NA
brfarmers$net12[which(!(brfarmers$net12 %in% surveyed))]  <- NA
brfarmers$net13[which(!(brfarmers$net13 %in% surveyed))]  <- NA

# Reshaping and droping lost
brfarmers.long <- reshape(
  brfarmers, v.names= "net",
  varying = c("net31", "net32", "net33", "net21", "net22", "net23", "net11", "net12", "net13"),
  timevar = "level", idvar="id", direction="long", drop = c("idold"))

# brfarmers.long <- subset(brfarmers.long, !is.na(net))

length(unique(unlist(subset(brfarmers.long, select=c(id,net)))))

library(netdiffuseR)

# Creating the graph object
graph <- with(brfarmers.long, edgelist_to_adjmat(cbind(id, net), undirected=TRUE, use.incomplete=FALSE, t=19))
used.vertex <- rownames(graph)

# Naming the array
dimnames(graph) <- list(used.vertex, used.vertex, 1948:1966)

# Average indegree
dg <- netdiffuseR::degree(graph, "indegree", undirected = FALSE)
dg <- rowMeans(dg + 1)
dg <- dg/max(dg)

# Difussion
toa <- brfarmers$yr[brfarmers$id %in% used.vertex]
adopt <- toa_mat(toa)
plot_diffnet(graph, adopt$cumadopt, displayisolates = FALSE, displaylabels=FALSE, mai = c(0,0,0,0), vertex.cex = 2)

# Infection
x <- plot_infectsuscep(graph, toa, K=10, logscale = TRUE, bins=20)

# Threshold
expo <- exposure(graph, adopt$cumadopt)
x <- plot_threshold(graph, expo, toa, undirected = FALSE, vertex.cex = 1/5)
x$fitted <- loess(threshold~jit, x, parametric = FALSE)$fitted
lines(fitted~jit, x[order(x$jit),], lwd=3, col="black")
