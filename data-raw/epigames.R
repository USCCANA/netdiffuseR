# data-raw/epigames.R
# Pre-processing script for the EpiGames Raw Dataset

rm(list = ls())

# The raw data consists of an attributes data frame and an hourly edgelist,
# both using consistent node IDs (1-594).
load("data-raw/epigames_hourly.rda")

epigames <- epigames_hourly

# Save compressed raw data
usethis::use_data(epigames, overwrite = TRUE, compress = "xz")
