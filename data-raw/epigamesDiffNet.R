# data-raw/epigamesDiffNet.R
# Generating the dynamic diffnet object using netdiffuseR + collapse_timeframes()

rm(list = ls())
library(netdiffuseR)

# Load the base raw dataset created in data-raw/epigames.R (hourly resolution)
load("data/epigames.rda")

attrs <- epigames$attributes
edges <- epigames$edgelist

# Collapse hourly edgelist (hours 0-338) into daily windows (days 1-15)
source("R/collapse_timeframes.R")

daily_edgelist <- collapse_timeframes(
  edgelist   = edges,
  ego        = "sender",
  alter      = "receiver",
  timevar    = "time",
  weightvar  = "weight",
  window_size = 24
)

# Build daily adjacency matrices
adjmat <- edgelist_to_adjmat(
  daily_edgelist[, c("sender", "receiver")],
  w  = daily_edgelist$weight,
  t0 = daily_edgelist$time,
  keep.isolates = TRUE,
  multiple      = TRUE
)

max_t <- max(daily_edgelist$time, na.rm = TRUE)

# Prepare TOA vector: real adoption times from attrs, NA for non-adopters
toa_vec <- stats::setNames(attrs$toa, as.character(attrs$id))

epigamesDiffNet <- as_diffnet(
  adjmat,
  toa                 = toa_vec,
  vertex.static.attrs = attrs,
  t0 = 1,
  t1 = max_t
)

usethis::use_data(epigamesDiffNet, overwrite = TRUE, compress = "xz")
