# data-raw/epigamesDiffNet.R
# Generating the dynamic diffnet object using netdiffuseR base methodology

rm(list = ls())
library(netdiffuseR)

# Load the base raw dataset created in data-raw/epigames.R
load("data/epigames_raw.rda")

attrs <- epigames_raw$attributes
edges <- epigames_raw$edgelist

# Creating dynamic graph via edgelist_to_adjmat (edges containing sender, receiver, time, weight)
adjmat <- edgelist_to_adjmat(
  edges[, c("sender", "receiver")],
  w = edges$weight,
  t = edges$time,
  keep.isolates = TRUE
)

# Coercing into a diffnet object
epigamesDiffNet <- as_diffnet(
  adjmat,
  toa = attrs$toa,
  vertex.attr = attrs
)

# Exporting formatted dataset for inclusion in the package
save(epigamesDiffNet, file = "data/epigamesDiffNet.rda", compress = "xz")

message("diffnet object successfully created and exported to data/epigamesDiffNet.rda")
