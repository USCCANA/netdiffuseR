rm(list=ls())
load("data/fakesurvey.rdata")

diffnet <- survey_to_diffnet(
  dat=fakesurvey,
  idvar = "id",
  netvars = c("net1", "net2", "net3"),
  toavar = "toa",
  groupvar="group",
  undirected = FALSE,
  multiple=TRUE)

fakeEdgelist <- data.frame(adjmat_to_edgelist(diffnet$graph,undirected = FALSE))
fakeEdgelist[,1] <- factor(diffnet$meta$ids[fakeEdgelist[,1]], diffnet$meta$ids)
fakeEdgelist[,2] <- factor(diffnet$meta$ids[fakeEdgelist[,2]], diffnet$meta$ids)
fakeEdgelist[,4] <- as.integer(fakeEdgelist[,4] )

# Creating a unique
fakeEdgelist$time <- NULL
fakeEdgelist <- unique(fakeEdgelist)

save(fakeEdgelist, file="data/fakeEdgelist.rdata")
