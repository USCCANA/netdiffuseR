# data-raw/epigames.R
# Pre-processing script for the Epi Games base dataset

rm(list = ls())

# Setting it up to use data/165/ as the source of raw parsed files for this repository mapping
attributes <- read.csv("data/165/epigames_attributes.csv", stringsAsFactors = FALSE)
edgelist <- read.csv("data/165/epigames_edgelist.csv", stringsAsFactors = FALSE)

# Package them into the expected raw struct
epigames_raw <- list(
  attributes = attributes,
  edgelist = edgelist
)

# Save as .rda compressed using xz for CRAN compliance
save(epigames_raw, file = "data/epigames_raw.rda", compress = "xz")

message("Data successfully compiled to data/epigames_raw.rda")
