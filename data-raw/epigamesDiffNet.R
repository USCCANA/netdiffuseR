# data-raw/epigamesDiffNet.R
# Generating the daily diffnet object from epigames using collapse_timeframes()
# Issue #75: Now includes vertex.dyn.attrs (mask, med, quarantine per day)
#
# Run after data-raw/epigames.R has built data/epigames.rda.

rm(list = ls())
library(netdiffuseR)

# ---------------------------------------------------------------------------
# 1. Load the base epigames dataset (with dynamic attrs)
# ---------------------------------------------------------------------------
load("data/epigames.rda")

attrs    <- epigames$attributes   # 594 x 6: id, toa, qyes_total, qno_total, mask_prop, med_prop
edges    <- epigames$edgelist     # hourly edgelist: sender, receiver, time (0-338), weight
dyn_long <- epigames$dyn_attrs    # long format: id, hour (0-338), mask, med, quarantine

# ---------------------------------------------------------------------------
# 2. Collapse hourly edgelist into 15 daily windows via collapse_timeframes()
# ---------------------------------------------------------------------------
WINDOW_SIZE <- 24   # hours per day
N_DAYS      <- 15

daily_edgelist <- collapse_timeframes(
  edgelist    = edges,
  ego         = "sender",
  alter       = "receiver",
  timevar     = "time",
  weightvar   = "weight",
  window_size = WINDOW_SIZE,
  binarize    = TRUE,
  cumulative  = TRUE,
  symmetric   = TRUE
)

cat("Daily edgelist: ", nrow(daily_edgelist), "rows, time range:", 
    range(daily_edgelist$time), "\n")

# Build adjacency matrices
adjmat <- edgelist_to_adjmat(
  daily_edgelist[, c("sender", "receiver")],
  w             = daily_edgelist$weight,
  t0            = daily_edgelist$time,
  keep.isolates = TRUE,
  multiple      = TRUE
)

# ---------------------------------------------------------------------------
# 3. Build vertex.dyn.attrs: one data.frame per day (15 total)
#    Each data.frame: 594 rows, columns: mask, med, quarantine (daily means)
# ---------------------------------------------------------------------------
# Map hourly data to day index (day d = hours [(d-1)*24 .. d*24-1])
dyn_long$day <- (dyn_long$hour %/% WINDOW_SIZE) + 1  # 1-based day
dyn_long$day <- pmin(dyn_long$day, N_DAYS)            # clamp hour 336-338 to day 15

vertex_dyn <- lapply(1:N_DAYS, function(d) {
  sub <- dyn_long[dyn_long$day == d, ]
  
  # Aggregate per node: mean within each 24-hour window
  # (proportion of hours in that day where behavior was active)
  agg <- aggregate(
    cbind(mask, med, quarantine) ~ id,
    data = sub,
    FUN  = mean
  )
  
  # Sort by id to match the node ordering in the diffnet object
  agg <- agg[order(agg$id), ]
  rownames(agg) <- NULL
  
  # Return only the behavior columns (not id — diffnet uses position)
  agg[, c("mask", "med", "quarantine")]
})

# Sanity check: each element should be 594 rows x 3 cols
stopifnot(all(sapply(vertex_dyn, nrow) == 594))
stopifnot(all(sapply(vertex_dyn, ncol) == 3))

cat("vertex.dyn.attrs built: ", N_DAYS, "data.frames of",
    nrow(vertex_dyn[[1]]), "rows x", ncol(vertex_dyn[[1]]), "cols\n")
cat("  Day 1  — mean mask usage:", round(mean(vertex_dyn[[1]]$mask), 3),
    "  mean quarantine:", round(mean(vertex_dyn[[1]]$quarantine), 3), "\n")
cat("  Day 15 — mean mask usage:", round(mean(vertex_dyn[[15]]$mask), 3),
    "  mean quarantine:", round(mean(vertex_dyn[[15]]$quarantine), 3), "\n")

# ---------------------------------------------------------------------------
# 4. Prepare TOA vector
# ---------------------------------------------------------------------------
toa_vec <- stats::setNames(attrs$toa, as.character(attrs$id))

# ---------------------------------------------------------------------------
# 5. Assemble diffnet object
# ---------------------------------------------------------------------------
epigamesDiffNet <- as_diffnet(
  adjmat,
  toa                 = toa_vec,
  vertex.static.attrs = attrs,
  vertex.dyn.attrs    = vertex_dyn,
  t0 = 1,
  t1 = N_DAYS
)

cat("\nepigamesDiffNet summary:\n")
print(epigamesDiffNet)

# ---------------------------------------------------------------------------
# 6. Quick validation: dynamic exposure vs static exposure
# ---------------------------------------------------------------------------
cat("\nValidating exposure() with dynamic mask attrs...\n")
expo_static  <- exposure(
  epigamesDiffNet,
  attrs = matrix(
    rep(epigamesDiffNet$vertex.static.attrs$mask_prop, N_DAYS),
    nrow = 594, ncol = N_DAYS
  )
)
expo_dynamic <- exposure(epigamesDiffNet, attrs = "mask")

cor_val <- cor(as.vector(expo_static), as.vector(expo_dynamic), use = "complete.obs")
cat("  Correlation static vs dynamic mask exposure:", round(cor_val, 4), "\n")
cat("  (Should be < 1.0, confirming dynamic attrs add new information)\n")

# ---------------------------------------------------------------------------
# 7. Save
# ---------------------------------------------------------------------------
usethis::use_data(epigamesDiffNet, overwrite = TRUE, compress = "xz")
cat("\nSaved: data/epigamesDiffNet.rda\n")
