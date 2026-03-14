#' Collapse Timeframes in a Longitudinal Edgelist
#'
#' @description
#' Allows users to take a high-resolution or continuous-time longitudinal
#' edgelist and dynamically collapse or discretize it into larger time windows.
#' The output is a shorter, aggregated edgelist ready to be passed into
#' \code{[edgelist_to_adjmat]} or \code{[as_diffnet]}.
#'
#' @param edgelist A \code{data.frame} representing the longitudinal edgelist.
#' @param ego Character scalar. Name of the column representing the ego (sender).
#' @param alter Character scalar. Name of the column representing the alter (receiver).
#' @param timevar Character scalar. Name of the column representing the time variable.
#' @param weightvar Character scalar or \code{NULL}. Name of the column representing
#'   the edge weight. If \code{NULL}, the function tallies the number of interactions
#'   within the time window as the weight.
#' @param window_size Numeric scalar. The size of the time window to collapse into.
#' @param time_format Character scalar or \code{NULL}. If the time variable is a
#'   character or factor, the format passed to \code{as.POSIXct}.
#'   For example, \code{"\%d-\%m-\%Y \%H:\%M"}.
#' @param relative_time Logical scalar. If \code{TRUE}, normalizes the binned
#'   times into a strict integer sequence starting at 1 (1, 2, 3...).
#'
#' @param binarize Logical scalar. If \code{TRUE}, sets all resulting edge weights to 1.
#' @param cumulative Logical scalar. If \code{TRUE}, edges from previous time windows
#'   are carried over to subsequent windows. 
#' @param symmetric Logical scalar. If \code{TRUE}, the resulting graph will be 
#'   symmetrized (i.e., if an edge A->B exists, an edge B->A is added).
#'
#' @return A \code{data.frame} with 4 columns: the ego, the alter, the new collapsed
#'   discrete time, and the aggregated weight.
#'
#' @export
#' @examples
#' \dontrun{
#' # Load the package's hourly dataset
#' load(system.file("data/epigames_raw.rda", package = "netdiffuseR"))
#' 
#' # Collapse the hourly edgelist into a daily edgelist (window_size = 24)
#' daily_edgelist <- collapse_timeframes(
#'   edgelist = epigames_raw$edgelist,
#'   timevar = "time",
#'   weightvar = "weight",
#'   window_size = 24
#' )
#' head(daily_edgelist)
#' }
collapse_timeframes <- function(
  edgelist, 
  ego = "sender", 
  alter = "receiver", 
  timevar = "time",
  weightvar = NULL,
  window_size = 1,
  time_format = NULL,
  relative_time = TRUE,
  binarize = FALSE,
  cumulative = FALSE,
  symmetric = FALSE
) {
  
  # Step 1: Time Column Parsing
  time_raw <- edgelist[[timevar]]
  
  if (is.character(time_raw) || is.factor(time_raw)) {
    if (!is.null(time_format)) {
      time_raw <- as.numeric(as.POSIXct(as.character(time_raw), format = time_format))
    } else {
      time_raw <- as.numeric(as.POSIXct(as.character(time_raw)))
    }
  } else if (!is.numeric(time_raw) && !is.integer(time_raw)) {
    # e.g., Date or POSIXct already
    time_raw <- as.numeric(time_raw)
  }
  
  # Check for NAs after conversion
  if (any(is.na(time_raw))) {
    warning("There are NA values in the parsed time variable.")
  }
  
  # Step 2: Binning / Window Creation
  t_min <- min(time_raw, na.rm = TRUE)
  # Adding a tiny offset so min time doesn't fall out of bounds or shift unnecessarily
  discrete_time <- ceiling((time_raw - t_min + 1e-9) / window_size)
  # Ensure the minimum index is 1 at this stage
  min_dt <- min(discrete_time, na.rm = TRUE)
  if (min_dt < 1) {
    discrete_time <- discrete_time - min_dt + 1
  }
  
  # Step 3: Handling relative_time
  if (relative_time) { # e.g. strict sequence 1, 2, 3
    sorted_unique_times <- sort(unique(discrete_time[!is.na(discrete_time)]))
    time_map <- stats::setNames(seq_along(sorted_unique_times), sorted_unique_times)
    discrete_time <- unname(time_map[as.character(discrete_time)])
  }
  
  # Create a working data frame to hold things
  work_df <- data.frame(
    ego_col = edgelist[[ego]],
    alter_col = edgelist[[alter]],
    time_col = discrete_time
  )
  
  # Step 4: Aggregation
  if (is.null(weightvar)) {
    work_df$weight_col <- 1
  } else {
    work_df$weight_col <- edgelist[[weightvar]]
  }
  
  # Remove rows with NAs in essential grouping variables
  work_df <- work_df[!is.na(work_df$ego_col) & !is.na(work_df$alter_col) & !is.na(work_df$time_col), ]
  
  agg_df <- stats::aggregate(
    weight_col ~ ego_col + alter_col + time_col, 
    data = work_df, 
    FUN = sum, 
    na.rm = TRUE
  )
  
  # Step 5: Output with 4 clean columns
  weight_col_name <- if (is.null(weightvar)) "weight" else weightvar
  colnames(agg_df) <- c(ego, alter, timevar, weight_col_name)
  
  # Step 6: Post-aggregation processing
  
  # 6.1 Binarize
  if (binarize) {
    agg_df[[weight_col_name]] <- 1
  }
  
  # 6.2 Symmetrize
  if (symmetric) {
    rev_df <- agg_df
    rev_df[[ego]] <- agg_df[[alter]]
    rev_df[[alter]] <- agg_df[[ego]]
    
    # Combine and de-duplicate (in case they already existed symmetrically)
    agg_df <- unique(rbind(agg_df, rev_df))
  }
  
  # 6.3 Cumulative
  if (cumulative) {
    all_periods <- sort(unique(agg_df[[timevar]]))
    if (length(all_periods) > 1) {
      cumulative_el <- agg_df[agg_df[[timevar]] == all_periods[1], ]
      for (t_idx in 2:length(all_periods)) {
        t <- all_periods[t_idx]
        current <- agg_df[agg_df[[timevar]] == t, ]
        prev <- cumulative_el[cumulative_el[[timevar]] == all_periods[t_idx - 1], ]
        if (nrow(prev) > 0) {
          prev[[timevar]] <- t
        }
        # Combine current window with previous accumulated edges and de-duplicate
        combined <- unique(rbind(current, prev))
        cumulative_el <- rbind(cumulative_el, combined)
      }
      agg_df <- cumulative_el
    }
  }
  
  # Apply standard sort for consistent outputs: time, ego, alter
  order_idx <- order(agg_df[[timevar]], agg_df[[ego]], agg_df[[alter]])
  agg_df <- agg_df[order_idx, ]
  rownames(agg_df) <- NULL
  
  return(agg_df)
}
