# data-raw/epigamesDiffNet.R
# Generating the daily diffnet object from epigames using collapse_timeframes()
# Run after data-raw/epigames.R has built data/epigames.rda.

rm(list = ls())
library(netdiffuseR)

load("data/epigames.rda")

attrs <- epigames$attributes # 594 x 6: id, toa, qyes_total, qno_total, mask_prop, med_prop
edges <- epigames$edgelist # hourly edgelist: sender, receiver, time (0-338), weight
dyn_long <- epigames$dyn_attrs # long format: id, hour (0-338), mask, med, quarantine

# Collapse hourly edgelist into 15 daily windows via collapse_timeframes()
WINDOW_SIZE <- 24
N_DAYS <- 15

dyn_long$day <- (dyn_long$hour %/% WINDOW_SIZE) + 1
dyn_long$day <- pmin(dyn_long$day, N_DAYS) # day mapping

daily_edgelist <- collapse_timeframes(
  edgelist    = edges,
  ego         = "sender",
  alter       = "receiver",
  timevar     = "time",
  weightvar   = "weight",
  window_size = WINDOW_SIZE,
  binarize    = FALSE,
  cumulative  = FALSE,
  symmetric   = TRUE
)

# Build adjacency matrices
adjmat <- edgelist_to_adjmat(
  daily_edgelist[, c("sender", "receiver")],
  w             = daily_edgelist$weight,
  t0            = daily_edgelist$time,
  t1            = daily_edgelist$time,
  keep.isolates = TRUE,
  multiple      = TRUE
)

# Build vertex.dyn.attrs: one data.frame per day (15 total)
# Each data.frame: 594 rows, columns: mask, med, quarantine

vertex_dyn <- lapply(1:N_DAYS, function(d) {
  sub <- dyn_long[dyn_long$day == d, ]

  # Aggregate per node: mean within each 24-hour window
  agg <- aggregate(
    cbind(mask, med, quarantine) ~ id,
    data = sub,
    FUN  = mean
  )

  # Sort by id to match the node ordering in the diffnet object
  agg <- agg[order(agg$id), ]
  rownames(agg) <- NULL

  # Return only the behavior columns
  agg[, c("mask", "med", "quarantine")]
})

# Prepare TOA vector
toa_vec <- stats::setNames(attrs$toa, as.character(attrs$id))

# Assemble diffnet object
epigamesDiffNet <- as_diffnet(
  adjmat,
  toa = toa_vec,
  vertex.static.attrs = attrs,
  vertex.dyn.attrs = vertex_dyn,
  t0 = 1,
  t1 = N_DAYS
)

# Save
usethis::use_data(epigamesDiffNet, overwrite = TRUE, compress = "xz")
