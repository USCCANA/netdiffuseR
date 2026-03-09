# data-raw/epigamesDiffNet.R
# Generating the dynamic diffnet object using netdiffuseR base methodology

rm(list = ls())
library(netdiffuseR)

# Load the base raw dataset created in data-raw/epigames.R (which is hourly)
load("data/epigames.rda")

attrs <- epigames$attributes
edges <- epigames$edgelist

# We need to collapse the hourly times into daily windows
# since diffnet objects often represent more macro periods
source("R/collapse_timeframes.R") 

daily_edgelist <- collapse_timeframes(
  edgelist = edges,
  ego = "sender",
  alter = "receiver",
  timevar = "time",
  weightvar = "weight",
  window_size = 24
)

# Creating dynamic graph via edgelist_to_adjmat
adjmat <- edgelist_to_adjmat(
  daily_edgelist[, c("sender", "receiver")],
  w = daily_edgelist$weight,
  t0 = daily_edgelist$time,
  keep.isolates = TRUE,
  multiple = TRUE
)

# Right censoring rule
# Let's check max time:
max_t <- max(daily_edgelist$time, na.rm = TRUE)

if (!"toa" %in% colnames(attrs)) {
  # If there is no TOA column, all nodes are non-adopters for this dataset.
  attrs$toa <- max_t + 1
} else if (any(is.na(attrs$toa))) {
  attrs$toa[is.na(attrs$toa)] <- max_t + 1
}

# Ensure attrs has matching ids with the edgelist.
# The graph has nodes numbered 1 to 594.
graph_nodes <- as.character(sort(as.numeric(rownames(adjmat[[1]]))))

if (nrow(attrs) != length(graph_nodes)) {
  attrs <- attrs[1:length(graph_nodes), , drop = FALSE]
  attrs$id <- as.integer(graph_nodes)
}
rownames(attrs) <- as.character(attrs$id)

# Create named toa vector using NA for non-adopters during diffnet creation
toa_vec <- stats::setNames(rep(NA, nrow(adjmat[[1]])), as.character(attrs$id))

# Let's verify and just format the diffnet cleanly
epigamesDiffNet <- as_diffnet(
  adjmat,
  toa = toa_vec, 
  vertex.static.attrs = attrs,
  t0 = 1,
  t1 = max_t
)

# Set non-adopters back to standard right-censoring rule after creation
epigamesDiffNet$toa[is.na(epigamesDiffNet$toa)] <- max_t + 1


# Exporting formatted dataset for inclusion in the package
usethis::use_data(epigamesDiffNet, overwrite = TRUE, compress = "xz")

message("diffnet object successfully created and exported to data/epigamesDiffNet.rda")