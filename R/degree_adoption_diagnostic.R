#' Degree and Time of Adoption Diagnostic
#'
#' Analyzes the correlation between in-degree, out-degree, and time of adoption
#' to identify whether opinion leaders were early adopters (supporters) or late
#' adopters (opposers).
#'
#' @param graph A `[diffnet()]` object or a graph data structure (classes include
#'   `array` (\eqn{n\times n \times T}{n*n*T}), `dgCMatrix` (sparse),
#'   `igraph`, etc.; see [netdiffuseR-graphs]).
#' @param toa Integer vector of length \eqn{n} (single behavior) or an \eqn{n\times Q}{n*Q}
#'   matrix (multi-behavior) with times of adoption. Required when `graph` is not a `diffnet`.
#' @param t0,t1 Optional integer scalars defining the first and last observed
#'   periods. If missing and `toa` is provided, `t0` defaults to 1
#'   and `t1` to `max(toa, na.rm=TRUE)`.
#' @param name Optional character scalars used only when coercing
#'   inputs into a `diffnet` object (passed to `new_diffnet`).
#' @param behavior Which behaviors to include when `toa` is a matrix (multi-diffusion).
#'   Can be `NULL` (all), a numeric index vector, or a character vector matching `colnames(toa)`.
#' @param combine Character scalar. How to combine multiple behaviors when `toa` is a matrix:
#'   - `"none"` (analyze each behavior separately)
#'   - `"pooled"` (stack rows across behaviors)
#'   - `"average"` (per-actor mean of TOA across selected behaviors)
#'   - `"earliest"` (per-actor minimum TOA)
#'   Ignored for single-behavior.
#' @param min_adopters Integer scalar. Minimum number of adopters required to compute correlations
#'   for any analysis cell (default 3).
#' @param degree_strategy Character scalar. How to aggregate degree measures across
#'   time periods:
#'   - `"mean"` (default): Average degree across all time periods
#'   - `"first"`: Degree in the first time period
#'   - `"last"`: Degree in the last time period
#' @param bootstrap Logical scalar. Whether to compute bootstrap confidence intervals.
#' @param R Integer scalar. Number of bootstrap replicates (default 1000).
#' @param conf.level Numeric scalar. Confidence level for bootstrap intervals (default 0.95).
#' @param valued Logical scalar. Whether to use edge weights in degree calculations.
#' @param ... Additional arguments passed on when coercing to `diffnet`.
#'
#' @details
#' This diagnostic function computes correlations between degree centrality measures
#' (in-degree and out-degree) and time of adoption. Positive correlations suggest
#' that central actors (opinion leaders) adopted early, while negative correlations
#' suggest they adopted late.
#'
#' When `bootstrap = TRUE`, the function uses the `boot` package to
#' compute bootstrap confidence intervals for the correlations.
#'
#' When `toa` is a matrix (multi-diffusion), degree vectors are computed once and
#' reused; the time of adoption is combined according to `combine`:
#' - `"none"`: computes separate results per behavior (see Value).
#' - `"pooled"`: stacks (actor, behavior) rows for adopters and runs a single analysis.
#' - `"average"`: one row per actor using the mean TOA of adopted behaviors.
#' - `"earliest"`: one row per actor using the minimum TOA of adopted behaviors.
#'
#' @return When analyzing a single behavior (or when `combine!="none"`), a list with:
#' \item{correlations}{Named numeric vector with correlations between in-degree/out-degree and time of adoption}
#' \item{bootstrap}{List with bootstrap results when `bootstrap = TRUE`, otherwise `NULL`}
#' \item{call}{The matched call}
#' \item{degree_strategy}{The degree aggregation strategy used}
#' \item{sample_size}{Number of rows included in the analysis (adopter rows)}
#' \item{combine}{`NULL` for single-behavior; otherwise the combination rule used.}
#'
#' When `combine="none"` with multiple behaviors, returns the same structure, except:
#' - `correlations` is a \eqn{2\times Q^*}{2 x Q*} matrix with rows `c("indegree_toa","outdegree_toa")`
#'   and one column per analyzed behavior.
#' - `bootstrap` is a named list with one entry per behavior (each like the single-behavior case), or `NULL` if `bootstrap=FALSE`.
#' - `sample_size` is an integer vector named by behavior.
#' - `combine` is `"none"`.
#'
#' @examples
#' # Basic usage with Korean Family Planning data
#' data(kfamilyDiffNet)
#' result_basics <- degree_adoption_diagnostic(kfamilyDiffNet, bootstrap = FALSE)
#' print(result_basics)
#'
#' # With bootstrap confidence intervals
#' result_boot <- degree_adoption_diagnostic(kfamilyDiffNet)
#' print(result_boot)
#'
#' # Different degree aggregation strategies
#' result_first <- degree_adoption_diagnostic(kfamilyDiffNet, degree_strategy = "first")
#' result_last  <- degree_adoption_diagnostic(kfamilyDiffNet, degree_strategy = "last")
#'
#' # Multi-diffusion (toy) ----------------------------------------------------
#' set.seed(999)
#' n <- 40; t <- 5; q <- 2
#' garr <- rgraph_ws(n, t, p=.3)
#' diffnet_multi <- rdiffnet(seed.graph = garr, t = t, seed.p.adopt = rep(list(0.1), q))
#'
#' # pooled (one combined analysis)
#' degree_adoption_diagnostic(diffnet_multi, combine = "pooled", bootstrap = FALSE)
#'
#' # per-behavior (matrix of correlations; one column per behavior)
#' degree_adoption_diagnostic(diffnet_multi, combine = "none", bootstrap = FALSE)
#'
#' @seealso `[dgr()]`, `[diffreg()]`, `[exposure()]`
#' @family statistics
#' @export
degree_adoption_diagnostic <- function(
  graph,
  degree_strategy = c("mean", "first", "last"),
  bootstrap = TRUE,
  R = 1000,
  conf.level = 0.95,
  toa = NULL,
  t0 = NULL, t1 = NULL,
  name = NULL,
  behavior = NULL,
  combine = c("none", "pooled", "average", "earliest"),
  min_adopters = 3,
  valued = getOption("diffnet.valued", FALSE),
  ...
) {
  # Check that bootstrap is a logical scalar
  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("'bootstrap' must be a logical scalar")
  }

  # Check that R is a positive integer
  if (!is.numeric(R) || length(R) != 1 || is.na(R) || R < 1 || R != as.integer(R)) {
    stop("'R' must be a positive integer")
  }

  # Check that conf.level is a numeric scalar between 0 and 1
  if (!is.numeric(conf.level) || length(conf.level) != 1 || is.na(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be between 0 and 1")
  }

  # Match and validate arguments
  degree_strategy <- match.arg(degree_strategy)
  combine <- match.arg(combine)

  # Store the original call
  original_call <- match.call()

  # Handle input processing and validation
  graph_and_toa <- process_graph_input(graph, toa, t0, t1, name, ...)
  graph <- graph_and_toa$graph
  toa <- graph_and_toa$toa

  # Get the number of behaviors
  if (is.matrix(toa)) {
    Q <- ncol(toa)
    behavior_names <- if (length(colnames(toa))) colnames(toa) else paste0("B", seq_len(Q))
  } else {
    Q <- 1
    behavior_names <- NULL
  }

  # Filter behaviors if requested
  if (Q > 1 && !is.null(behavior)) {
    if (is.character(behavior)) {
      if (is.null(colnames(toa))) {
        stop("Cannot use character behavior selection without column names in toa")
      }
      behavior_indices <- match(behavior, colnames(toa))
      if (any(is.na(behavior_indices))) {
        stop("Some behavior names not found in colnames(toa): ",
             paste(behavior[is.na(behavior_indices)], collapse = ", "))
      }
    } else if (is.numeric(behavior)) {
      behavior_indices <- behavior
      if (any(behavior_indices < 1 | behavior_indices > Q)) {
        stop("behavior index out of range: must be between 1 and ", Q)
      }
    } else {
      stop("behavior must be NULL, numeric indices, or character names")
    }

    toa <- toa[, behavior_indices, drop = FALSE]
    Q <- ncol(toa)
    behavior_names <- if (length(colnames(toa))) colnames(toa) else paste0("B", seq_len(Q))
  }

  # Compute degree measures
  degrees <- compute_degree_measures(graph, degree_strategy, valued)

  # Handle multi-diffusion analysis: for combine = "none", always analyze per behavior and do not error for too few adopters
  if (combine == "none" && is.matrix(toa)) {
    return(analyze_multi_behaviors_separately(
      degrees, toa, min_adopters, bootstrap, R, conf.level,
      behavior_names, degree_strategy, original_call, graph
    ))
  }

  # Single behavior or combined analysis
  combined_data <- prepare_combined_data(degrees, toa, combine, min_adopters, Q)

  if (nrow(combined_data) < min_adopters) {
    stop("Insufficient adopters for correlation analysis. (n=", nrow(combined_data),
      ", minimum = ", min_adopters, ").")
  }

  # Compute correlations
  correlations <- compute_correlations(combined_data)

  # Bootstrap if requested
  bootstrap_results <- if (bootstrap) {
    compute_bootstrap_results(combined_data, R, conf.level)
  } else {
    NULL
  }

  # Determine if undirected (graph is always a diffnet here)
  undirected <- isTRUE(is_undirected(graph))

  # Return results
  structure(list(
    correlations = correlations,
    bootstrap = bootstrap_results,
    call = original_call,
    degree_strategy = degree_strategy,
    sample_size = nrow(combined_data),
    combine = if (Q > 1) combine else NULL,
    undirected = undirected
  ), class = "degree_adoption_diagnostic")
}

# Helper functions --------------------------------------------------------

process_graph_input <- function(graph, toa, t0, t1, name, ...) {
  if (inherits(graph, "diffnet")) {
    if (!is.null(toa)) {
      warning("toa argument ignored when graph is a diffnet object")
    }
    return(list(graph = graph, toa = graph$toa))
  }

  if (is.null(toa)) {
    stop("toa argument is required when graph is not a diffnet object")
  }

  # If graph is a list, ensure all elements are dgCMatrix
  if (is.list(graph)) {
    graph <- lapply(graph, function(g) {
      if (inherits(g, "dgCMatrix")) return(g)
      if (is.matrix(g)) return(as(Matrix::Matrix(g, sparse = TRUE), "dgCMatrix"))
      stop("All elements of the graph list must be matrices or dgCMatrix.")
    })
  }

  # If graph is a single static adjacency (matrix/dgCMatrix), expand it to T layers to avoid recycling warnings in new_diffnet
  if (inherits(graph, "dgCMatrix") || is.matrix(graph)) {
    if (is.null(t0)) t0 <- 1
    if (is.null(t1)) t1 <- max(toa, na.rm = TRUE)
    Tlen <- as.integer(t1 - t0 + 1L)
    if (is.finite(Tlen) && Tlen > 1L) {
      gmat <- if (inherits(graph, "dgCMatrix")) graph else as(Matrix::Matrix(graph, sparse = TRUE), "dgCMatrix")
      graph <- replicate(Tlen, gmat, simplify = FALSE)
    }
  }

  # Detect undirectedness on the raw input before building diffnet
  undirected_flag <- check_undirected_graph(graph)

  # Convert to diffnet
  if (is.null(t0)) t0 <- 1
  if (is.null(t1)) t1 <- max(toa, na.rm = TRUE)

  graph <- new_diffnet(
    graph, toa,
    t0 = t0, t1 = t1,
    name = name,
    undirected = undirected_flag,
    ...
  )
  return(list(graph = graph, toa = toa))
}

compute_degree_measures <- function(graph, degree_strategy, valued) {
  if (degree_strategy == "mean") {
    indegree <- rowMeans(dgr(graph, cmode = "indegree", valued = valued), na.rm = TRUE)
    outdegree <- rowMeans(dgr(graph, cmode = "outdegree", valued = valued), na.rm = TRUE)
  } else {
    deg_matrix <- dgr(graph, valued = valued)
    if (length(dim(deg_matrix)) == 3) {
      # Dynamic case
      if (degree_strategy == "first") {
        indegree <- deg_matrix[, 1, "indegree"]
        outdegree <- deg_matrix[, 1, "outdegree"]
      } else if (degree_strategy == "last") {
        last_time <- dim(deg_matrix)[2]
        indegree <- deg_matrix[, last_time, "indegree"]
        outdegree <- deg_matrix[, last_time, "outdegree"]
      }
    } else if (length(dim(deg_matrix)) == 2) {
      # Static case: check for column names, else use position
      cn <- colnames(deg_matrix)
      if (!is.null(cn) && all(c("indegree", "outdegree") %in% cn)) {
        indegree <- deg_matrix[, "indegree"]
        outdegree <- deg_matrix[, "outdegree"]
      } else if (ncol(deg_matrix) >= 2) {
        indegree <- deg_matrix[, 1]
        outdegree <- deg_matrix[, 2]
      } else {
        stop("Degree matrix does not have expected columns for static graph.")
      }
    } else {
      stop("Unexpected degree matrix dimensions in compute_degree_measures.")
    }
  }

  list(indegree = indegree, outdegree = outdegree)
}

analyze_multi_behaviors_separately <- function(degrees, toa, min_adopters, bootstrap, R, conf.level, behavior_names, degree_strategy, original_call, graph) {
  Q <- ncol(toa)

  # Initialize results containers
  correlations_matrix <- matrix(NA_real_, nrow = 2, ncol = Q)
  rownames(correlations_matrix) <- c("indegree_toa", "outdegree_toa")
  colnames(correlations_matrix) <- behavior_names

  sample_sizes <- integer(Q)
  names(sample_sizes) <- behavior_names

  bootstrap_results <- if (bootstrap) vector("list", Q) else NULL
  if (bootstrap) names(bootstrap_results) <- behavior_names

  # Analyze each behavior
  for (q in seq_len(Q)) {
    toa_q <- toa[, q]
    adopters_q <- which(!is.na(toa_q))

    if (length(adopters_q) >= min_adopters) {
      data_q <- data.frame(
        indegree = degrees$indegree[adopters_q],
        outdegree = degrees$outdegree[adopters_q],
        toa = toa_q[adopters_q]
      )

      correlations_matrix[1, q] <- cor_safe(data_q$indegree, data_q$toa )
      correlations_matrix[2, q] <- cor_safe(data_q$outdegree, data_q$toa )
      sample_sizes[q] <- nrow(data_q)

      if (bootstrap) {
        bootstrap_results[[q]] <- compute_bootstrap_results(data_q, R, conf.level)
      }
    } else {
      sample_sizes[q] <- length(adopters_q)
    }
  }

  # Determine if undirected
  undirected <- if (inherits(graph, "diffnet")) {
    is_undirected(graph)
  } else {
    check_undirected_graph(graph)
  }

  structure(list(
    correlations = correlations_matrix,
    bootstrap = bootstrap_results,
    call = original_call,
    degree_strategy = degree_strategy,
    sample_size = sample_sizes,
    combine = "none",
    undirected = undirected
  ), class = "degree_adoption_diagnostic")
}

prepare_combined_data <- function(degrees, toa, combine, min_adopters, Q) {
  if (Q == 1 || combine == "pooled") {
    if (Q == 1) {
      adopters <- which(!is.na(toa))
      data.frame(
        indegree = degrees$indegree[adopters],
        outdegree = degrees$outdegree[adopters],
        toa = toa[adopters]
      )
    } else {
      # Pooled: stack all (actor, behavior) rows
      adopter_rows <- which(!is.na(as.vector(toa)))
      actor_indices <- ((adopter_rows - 1) %% nrow(toa)) + 1

      data.frame(
        indegree = degrees$indegree[actor_indices],
        outdegree = degrees$outdegree[actor_indices],
        toa = as.vector(toa)[adopter_rows]
      )
    }
  } else if (combine == "average") {
    # Average TOA across behaviors per actor
    toa_avg <- rowMeans(toa, na.rm = TRUE)
    adopters <- which(!is.na(toa_avg))

    data.frame(
      indegree = degrees$indegree[adopters],
      outdegree = degrees$outdegree[adopters],
      toa = toa_avg[adopters]
    )
  } else if (combine == "earliest") {
    # Earliest TOA across behaviors per actor
    toa_min <- apply(toa, 1, function(row) {
      if (all(is.na(row))) return(NA_real_)
      min(row, na.rm = TRUE)
    })
    toa_min[is.infinite(toa_min)] <- NA
    adopters <- which(!is.na(toa_min))

    data.frame(
      indegree = degrees$indegree[adopters],
      outdegree = degrees$outdegree[adopters],
      toa = toa_min[adopters]
    )
  }
}

compute_correlations <- function(data) {
  c(
    indegree_toa = cor_safe(data$indegree, data$toa),
    outdegree_toa = cor_safe(data$outdegree, data$toa)
  )
}

compute_bootstrap_results <- function(data, R, conf.level) {
  safe_bootstrap <- function(data, indices) {
    d <- data[indices, ]
    tryCatch({
      c(
        indegree_toa = cor_safe(d$indegree, d$toa ),
        outdegree_toa = cor_safe(d$outdegree, d$toa )
      )
    }, error = function(e) c(indegree_toa = NA, outdegree_toa = NA))
  }

  boot_result <- boot::boot(data, safe_bootstrap, R = R)

  # Observed correlations
  obs_indegree <- cor_safe(data$indegree, data$toa)
  obs_outdegree <- cor_safe(data$outdegree, data$toa)

  # Bias and standard error
  bias_indegree <- mean(boot_result$t[,1], na.rm = TRUE) - obs_indegree
  std_error_indegree <- sd(boot_result$t[,1], na.rm = TRUE)
  bias_outdegree <- mean(boot_result$t[,2], na.rm = TRUE) - obs_outdegree
  std_error_outdegree <- sd(boot_result$t[,2], na.rm = TRUE)

  # Confidence intervals
  indegree_ci <- tryCatch({
    ci <- boot::boot.ci(boot_result, conf = conf.level, type = "perc", index = 1)
    ci$percent[4:5]
  }, error = function(e) c(NA, NA))

  outdegree_ci <- tryCatch({
    ci <- boot::boot.ci(boot_result, conf = conf.level, type = "perc", index = 2)
    ci$percent[4:5]
  }, error = function(e) c(NA, NA))

  list(
    indegree = list(
      correlation = obs_indegree,
      bias = bias_indegree,
      std_error = std_error_indegree,
      conf_int = indegree_ci,
      conf_level = conf.level
    ),
    outdegree = list(
      correlation = obs_outdegree,
      bias = bias_outdegree,
      std_error = std_error_outdegree,
      conf_int = outdegree_ci,
      conf_level = conf.level
    ),
    R = R,
    boot_object = boot_result
  )
}

create_empty_result <- function(degree_strategy, original_call, combine, sample_size) {
  structure(list(
    correlations = c(indegree_toa = NA_real_, outdegree_toa = NA_real_),
    bootstrap = NULL,
    call = original_call,
    degree_strategy = degree_strategy,
    sample_size = sample_size,
    combine = combine,
    undirected = FALSE
  ), class = "degree_adoption_diagnostic")
}

check_undirected_graph <- function(graph) {
  if (is.list(graph)) {
    return(all(sapply(graph, function(g) isSymmetric(as.matrix(g)))))
  }
  if (is.array(graph) && length(dim(graph)) == 3) {
    return(all(sapply(seq_len(dim(graph)[3]), function(t) isSymmetric(as.matrix(graph[,,t])))))
  }
  if (is.matrix(graph)) {
    return(isSymmetric(as.matrix(graph)))
  }
  FALSE
}

# Print method ------------------------------------------------------------

#' @export
print.degree_adoption_diagnostic <- function(x, ...) {
  cat("Degree and Time of Adoption Diagnostic\n")
  cat("======================================\n\n")

  # Basic info
  cat(sprintf("Degree aggregation strategy: %s \n", x$degree_strategy))

  # Sample size info
  if (is.null(names(x$sample_size))) {
    cat(sprintf("Sample size (adopters only): %d \n", x$sample_size))
  } else {
    cat("Sample sizes (adopters only):\n")
    beh_names <- names(x$sample_size)
    for (j in seq_along(x$sample_size)) {
      cat(sprintf("    - %s: %d\n", if (length(beh_names)) beh_names[j] else paste0("B", j), x$sample_size[j]))
    }
    cat("\n")
  }

  undirected <- isTRUE(x$undirected)

  # Single behavior or combined analysis
  if (is.vector(x$correlations)) {
    print_single_behavior_results(x, undirected)
  } else {
    print_multi_behavior_results(x, undirected)
  }

  invisible(x)
}

print_single_behavior_results <- function(x, undirected) {
  # Extract correlation values
  indeg_r <- x$correlations[["indegree_toa"]]
  outdeg_r <- x$correlations[["outdegree_toa"]]

  # Print correlations
  cat("\nCorrelations:\n")
  if (undirected) {
    deg_r <- indeg_r  # For undirected graphs, in-degree = out-degree = degree
    cat(sprintf("  Degree - Time of Adoption: %.3f\n", deg_r))
  } else {
    cat(sprintf("  In-degree  - Time of Adoption: %.3f\n", indeg_r))
    cat(sprintf("  Out-degree - Time of Adoption: %.3f\n", outdeg_r))
  }

  # Interpretation
  cat("\nInterpretation:\n")

  if (!is.null(x$bootstrap)) {
    bootstrap_data <- x$bootstrap
    deg_ci <- if (undirected && !is.null(bootstrap_data$indegree$conf_int)) {
      bootstrap_data$indegree$conf_int
    } else NULL
    indeg_ci <- if (!is.null(bootstrap_data$indegree$conf_int)) {
      bootstrap_data$indegree$conf_int
    } else NULL
    outdeg_ci <- if (!is.null(bootstrap_data$outdegree$conf_int)) {
      bootstrap_data$outdegree$conf_int
    } else NULL
    lvl <- if (!is.null(bootstrap_data$indegree$conf_level)) {
      bootstrap_data$indegree$conf_level * 100
    } else NA_real_

    if (undirected) {
      explain_degree_correlation("Degree", deg_r, deg_ci, lvl_arg = lvl)
    } else {
      explain_degree_correlation("In-degree", indeg_r, indeg_ci, lvl_arg = lvl)
      explain_degree_correlation("Out-degree", outdeg_r, outdeg_ci, lvl_arg = lvl)
    }
  } else {
    if (undirected) {
      explain_degree_correlation("Degree", deg_r, NULL)
    } else {
      explain_degree_correlation("In-degree", indeg_r, NULL)
      explain_degree_correlation("Out-degree", outdeg_r, NULL)
    }
  }
}

print_multi_behavior_results <- function(x, undirected) {
  correlations_matrix <- x$correlations
  Q <- ncol(correlations_matrix)
  beh_names <- colnames(correlations_matrix)

  # Print correlations matrix
  cat("\nCorrelations:\n")
  if (undirected) {
    cat("  Degree - Time of Adoption:\n")
    deg_row <- correlations_matrix["indegree_toa", ]
    for (j in seq_len(Q)) {
      bname <- if (length(beh_names)) beh_names[j] else paste0("B", j)
      cat(sprintf("    [%s]: %.3f\n", bname, deg_row[j]))
    }
  } else {
    cat("  In-degree - Time of Adoption:\n")
    indeg_row <- correlations_matrix["indegree_toa", ]
    for (j in seq_len(Q)) {
      bname <- if (length(beh_names)) beh_names[j] else paste0("B", j)
      cat(sprintf("    [%s]: %.3f\n", bname, indeg_row[j]))
    }
    cat("  Out-degree - Time of Adoption:\n")
    outdeg_row <- correlations_matrix["outdegree_toa", ]
    for (j in seq_len(Q)) {
      bname <- if (length(beh_names)) beh_names[j] else paste0("B", j)
      cat(sprintf("    [%s]: %.3f\n", bname, outdeg_row[j]))
    }
  }

  # Interpretation per behavior
  cat("\nInterpretation:\n")
  for (j in seq_len(Q)) {
    bname <- if (length(beh_names)) beh_names[j] else paste0("B", j)
    r_in <- correlations_matrix["indegree_toa", j]
    r_out <- correlations_matrix["outdegree_toa", j]

    bootstrap_data <- if (!is.null(x$bootstrap)) x$bootstrap[[j]] else NULL
    deg_ci <- if (undirected && !is.null(bootstrap_data) && !is.null(bootstrap_data$indegree$conf_int)) {
      bootstrap_data$indegree$conf_int
    } else NULL
    indeg_ci <- if (!is.null(bootstrap_data) && !is.null(bootstrap_data$indegree$conf_int)) {
      bootstrap_data$indegree$conf_int
    } else NULL
    outdeg_ci <- if (!is.null(bootstrap_data) && !is.null(bootstrap_data$outdegree$conf_int)) {
      bootstrap_data$outdegree$conf_int
    } else NULL
    lvl <- if (!is.null(bootstrap_data) && !is.null(bootstrap_data$indegree$conf_level)) {
      bootstrap_data$indegree$conf_level * 100
    } else NA_real_

    cat(sprintf("  [%s]\n", bname))
    if (undirected) {
      explain_degree_correlation("Degree", r_in, deg_ci, lvl_arg = lvl)
    } else {
      explain_degree_correlation("In-degree", r_in, indeg_ci, lvl_arg = lvl)
      explain_degree_correlation("Out-degree", r_out, outdeg_ci, lvl_arg = lvl)
    }
  }
}

# Helper function for explaining correlations
explain_degree_correlation <- function(label, r, ci, lvl_arg = NA_real_, thr = 0.10) {
  # Handle NA correlations gracefully
  if (is.na(r)) {
    cat(sprintf(
      "  %s: Weak relationship between centrality and adoption timing:\n             r is NA; no CI.\n",
      label
    ))
    return(invisible())
  }

  abs_big <- abs(r) > thr
  degree_term <- switch(label,
    "In-degree" = "in-degree",
    "Out-degree" = "out-degree",
    "degree"
  )

  if (is.null(ci)) {
    format_interpretation_no_ci(label, r, abs_big, degree_term, thr)
  } else {
    format_interpretation_with_ci(label, r, ci, abs_big, degree_term, thr, lvl_arg)
  }
}

format_interpretation_no_ci <- function(label, r, abs_big, degree_term, thr) {
  if (!abs_big) {
    cat(sprintf("  %s: Weak relationship between %s and adoption timing:\n             |r| \\u2264 %.1f; no CI.\n",
                label, degree_term, thr))
  } else if (r > 0) {
    cat(sprintf("  %s: Central actors (high %s) tended to adopt early (supporters):\n             |r| > %.1f; no CI.\n",
                label, degree_term, thr))
  } else {
    cat(sprintf("  %s: Central actors (high %s) tended to adopt late (opposers):\n             |r| > %.1f; no CI.\n",
                label, degree_term, thr))
  }
}

format_interpretation_with_ci <- function(label, r, ci, abs_big, degree_term, thr, lvl_arg) {
  lvl_local <- if (!is.na(lvl_arg)) lvl_arg else 95
  ci_includes_zero <- (length(ci) >= 2) && is.finite(ci[1]) && is.finite(ci[2]) && (ci[1] <= 0 && ci[2] >= 0)

  if (!abs_big) {
    cat(sprintf("  %s: Weak relationship between %s and adoption timing; %s statistically supported:\n             |r| \\u2264 %.1f; CI (%.1f%%) %s 0.\n",
                label, degree_term,
                if (ci_includes_zero) "NOT" else "",
                thr, lvl_local,
                if (ci_includes_zero) "includes" else "excludes"))
  } else if (r > 0) {
    cat(sprintf("  %s: Central actors (high %s) tended to adopt early (supporters); %s statistically supported:\n             |r| > %.1f; CI (%.1f%%) %s 0.\n",
                label, degree_term,
                if (ci_includes_zero) "NOT" else "",
                thr, lvl_local,
                if (ci_includes_zero) "includes" else "excludes"))
  } else {
    cat(sprintf("  %s: Central actors (high %s) tended to adopt late (opposers); %s statistically supported:\n             |r| > %.1f; CI (%.1f%%) %s 0.\n",
                label, degree_term,
                if (ci_includes_zero) "NOT" else "",
                thr, lvl_local,
                if (ci_includes_zero) "includes" else "excludes"))
  }
}

# Safe correlation: returns NA (no warnings) if zero-variance or too few pairs
cor_safe <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}
