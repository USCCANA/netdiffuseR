rm(list=ls())
load("data/fakesurveyDyn.rdata")

diffnet <- survey_to_diffnet(
  dat=fakesurveyDyn,
  idvar = "id",
  netvars = c("net1", "net2", "net3"),
  toavar = "toa",
  groupvar="group",
  timevar = "time",
  undirected = FALSE,
  multiple=TRUE)

fakeDynEdgelist <- data.frame(adjmat_to_edgelist(diffnet$graph,undirected = FALSE))
fakeDynEdgelist[,1] <- factor(diffnet$meta$ids[fakeDynEdgelist[,1]], diffnet$meta$ids)
fakeDynEdgelist[,2] <- factor(diffnet$meta$ids[fakeDynEdgelist[,2]], diffnet$meta$ids)
fakeDynEdgelist[,4] <- as.integer(fakeDynEdgelist[,4] )

save(fakeDynEdgelist, file="data/fakeDynEdgelist.rdata")
