# data-raw/epigames.R
# Pre-processing script for the EpiGames Raw Dataset

rm(list = ls())

# The raw data consists of an attributes data frame and an hourly edgelist,
# both using consistent node IDs (1-594).
load("data-raw/epigames_hourly.rda")

# Load the hourly dynamic behavioral attributes
dyn_attrs_path <- "playground/epigames-stuff/epigames-analysis-copy/dynamic_attrs_hourly.csv"

dyn_attrs_hourly <- read.csv(dyn_attrs_path, stringsAsFactors = FALSE)

# Sanity checks
stopifnot(ncol(dyn_attrs_hourly) == 5) # id, hour, mask, med, quarantine
stopifnot(nrow(dyn_attrs_hourly) == 594 * 339) # 201,366 rows
stopifnot(all(dyn_attrs_hourly$id %in% 1:594))
stopifnot(all(dyn_attrs_hourly$hour %in% 0:338))

# Bundle into the epigames list (3 elements)
epigames <- list(
  attributes = epigames_hourly$attributes, # static, 594 x 6
  edgelist   = epigames_hourly$edgelist, # hourly, ~39k rows
  dyn_attrs  = dyn_attrs_hourly # dynamical attributes (long format)
)

# Save compressed .rda
usethis::use_data(epigames, overwrite = TRUE, compress = "xz")
