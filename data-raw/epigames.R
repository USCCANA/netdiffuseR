# data-raw/epigames.R
# Pre-processing script for the EpiGames Raw Dataset

rm(list=ls())

# The raw data is originally packaged from hourly resolution arrays.
# It consists of an attributes dataframe and an edgelist dataframe.
load("data-raw/epigames_hourly.rda")

# The data in the rda was saved as `epigames_hourly`. 
# We simply rename it to the package standard `epigames`
epigames <- epigames_hourly

# Save compressed raw data
usethis::use_data(epigames, overwrite = TRUE, compress = "xz")

message("Data successfully compiled to data/epigames.rda")
