#' Source-attribution kernels for \code{rdiffnet}'s lineage tracking
#'
#' A family of pluggable rules that decide, for each freshly adopted node
#' in an \code{\link{rdiffnet}} simulation, \emph{which} of its already
#' adopted neighbours infected it. Pass any of these as the
#' \code{source_attribution} argument of \code{\link{rdiffnet}}, or write
#' your own function that follows the same contract.
#'
#' @param target Integer scalar. Index of the freshly adopted node.
#' @param adopted_neighbours Integer vector of indices of neighbours of
#'   \code{target} that were already adopted at the previous time step.
#'   Pre-sorted by ascending time-of-adoption (earliest infector first).
#'   Empty when \code{target} is a seed.
#' @param weights Numeric vector of edge weights aligned with
#'   \code{adopted_neighbours} (same length, same order). Carries the
#'   non-zero entries of \code{sgraph[[time]][target, ]}. \code{NULL} or
#'   constant when the graph is unweighted.
#' @param time Integer scalar. Current simulation time step.
#' @param pars Named list. Passed verbatim from \code{adoption_pars} so
#'   user-defined attributors can read whatever they need.
#'
#' @return A single integer index — one of \code{adopted_neighbours},
#'   identifying the attributed source. \code{NA_integer_} when
#'   \code{adopted_neighbours} is empty (seed or spontaneously adopted
#'   node).
#'
#' @details
#' The contract is intentionally minimal. \code{rdiffnet()} pre-computes
#' \code{adopted_neighbours} and \code{weights} once per fresh adoption
#' and hands them sorted by toa, so user-defined attributors don't have
#' to query the simulation state themselves. The three kernels shipped
#' with the package:
#'
#' \describe{
#'  \item{\code{source_attribution_uniform}}{Samples uniformly across
#'    \code{adopted_neighbours}. The default when the user passes
#'    \code{TRUE} (deferred) and the simplest reasonable choice when
#'    nothing distinguishes the candidates.}
#'  \item{\code{source_attribution_weighted}}{Samples with probability
#'    proportional to \code{weights}. Falls back to uniform when the
#'    graph carries no weights (all entries equal). Appropriate for
#'    contact-network simulations where edge weights encode contact
#'    intensity.}
#'  \item{\code{source_attribution_earliest}}{Returns the
#'    earliest-infected adopted neighbour. Mirrors the heuristic that
#'    the playground's \code{derive_tree} used post-hoc; useful as a
#'    deterministic baseline.}
#' }
#'
#' @references
#' Lloyd-Smith, J. O., Schreiber, S. J., Kopp, P. E., & Getz, W. M. (2005).
#' Superspreading and the effect of individual variation on disease emergence.
#' \emph{Nature} 438:355-359. \doi{10.1038/nature04153}
#'
#' @examples
#' set.seed(2026)
#'
#' # Use a kernel directly inside rdiffnet():
#' dn <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
#'                seed.p.adopt = 0.1, stop.no.diff = FALSE,
#'                source_attribution = source_attribution_weighted)
#'
#' is.diffnet_epi(dn)              # TRUE — auto-promoted
#' nrow(transmission_tree(dn))      # one row per fresh adoption + seeds
#'
#' @name source_attribution
#' @author Aníbal Olivera M.
NULL

#' @rdname source_attribution
#' @export
source_attribution_uniform <- function(target, adopted_neighbours, weights,
                                       time, pars) {
  if (!length(adopted_neighbours)) return(NA_integer_)
  if (length(adopted_neighbours) == 1L) return(as.integer(adopted_neighbours[1L]))
  as.integer(sample(adopted_neighbours, size = 1L))
}

#' @rdname source_attribution
#' @export
source_attribution_weighted <- function(target, adopted_neighbours, weights,
                                        time, pars) {
  if (!length(adopted_neighbours)) return(NA_integer_)
  if (length(adopted_neighbours) == 1L) return(as.integer(adopted_neighbours[1L]))

  # No discriminating weights -> fall back to uniform silently. Users who
  # really want weighted attribution on an unweighted graph need to encode
  # the weights they care about into the simulation's graph.
  if (is.null(weights) || all(is.na(weights)) ||
      length(unique(weights)) == 1L) {
    return(as.integer(sample(adopted_neighbours, size = 1L)))
  }

  # Guard against zero / negative weights (sample() with prob requires non-
  # negative non-zero sums).
  w <- pmax(as.numeric(weights), 0)
  if (sum(w) == 0)
    return(as.integer(sample(adopted_neighbours, size = 1L)))

  as.integer(sample(adopted_neighbours, size = 1L, prob = w))
}

#' @rdname source_attribution
#' @export
source_attribution_earliest <- function(target, adopted_neighbours, weights,
                                        time, pars) {
  if (!length(adopted_neighbours)) return(NA_integer_)
  # By contract -adopted_neighbours- is sorted by ascending toa, so the
  # earliest infector is at position 1.
  as.integer(adopted_neighbours[1L])
}

# Internal: validate -source_attribution- and normalize it to a length-Q list.
# Returns a list of length num_of_behaviors with either NULL or a function per
# behaviour.
rdiffnet_normalize_source_attribution <- function(source_attribution,
                                                  num_of_behaviors) {
  if (is.null(source_attribution)) {
    return(rep(list(NULL), num_of_behaviors))
  }

  if (is.function(source_attribution)) {
    # Broadcast: same attributor for every behaviour.
    return(rep(list(source_attribution), num_of_behaviors))
  }

  if (is.list(source_attribution)) {
    if (length(source_attribution) != num_of_behaviors)
      stop("-source_attribution- list must have length ", num_of_behaviors,
           " (one entry per behaviour), got ", length(source_attribution), ".")
    bad <- which(!vapply(source_attribution,
                         function(z) is.null(z) || is.function(z),
                         logical(1)))
    if (length(bad))
      stop("Every -source_attribution[[q]]- must be NULL or a function. ",
           "Offending position(s): ", paste(bad, collapse = ", "), ".")
    return(source_attribution)
  }

  stop("-source_attribution- must be NULL, a function, or a length-Q list of ",
       "functions / NULLs.")
}

# Internal: assemble accumulated tree rows (list of one-row lists) into the
# canonical six-column transmission data.frame. -as_transmission_tree- /
# -as_diffnet_epi- will validate the result.
rdiffnet_tree_rows_to_df <- function(rows) {
  if (!length(rows)) {
    return(data.frame(
      date                 = integer(0),
      source               = integer(0),
      target               = integer(0),
      source_exposure_date = integer(0),
      virus_id             = integer(0),
      virus                = character(0),
      stringsAsFactors     = FALSE
    ))
  }
  data.frame(
    date                 = vapply(rows, `[[`, integer(1L),   "date"),
    source               = vapply(rows, function(r) as.integer(r$source),
                                  integer(1L)),
    target               = vapply(rows, `[[`, integer(1L),   "target"),
    source_exposure_date = vapply(rows,
                                  function(r) as.integer(r$source_exposure_date),
                                  integer(1L)),
    virus_id             = vapply(rows, `[[`, integer(1L),   "virus_id"),
    virus                = vapply(rows, `[[`, character(1L), "virus"),
    stringsAsFactors     = FALSE
  )
}
