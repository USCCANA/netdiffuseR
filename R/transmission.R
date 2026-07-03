# Transmission tree handling for diffnet objects.
#
# The $transmission slot is a list with the following elements:
#   - tree: data.frame with columns
#       date                  integer, period when the transmission happened
#       source                integer, row index of the infector in x (NA for seeds)
#       target                integer, row index of the infectee in x
#       source_exposure_date  integer, period when `source` was infected (NA for seeds)
#       virus_id              integer, optional virus identifier
#       virus                 character, optional virus label
#   - pars: list, free-form parameters/metadata associated with the tree.
# Each row represents one infection event (an edge in the transmission tree);
# the set of (source -> target) pairs forms the directed forest from which
# offspring distributions (Lloyd-Smith et al., 2005) and likelihood-based
# reproduction-number estimates (White & Pagano, 2008) are derived.

.transmission_cols <- c(
  "date", "source", "target", "source_exposure_date", "virus_id", "virus"
)

.empty_transmission_tree <- function() {
  data.frame(
    date                 = integer(0),
    source               = integer(0),
    target               = integer(0),
    source_exposure_date = integer(0),
    virus_id             = integer(0),
    virus                = character(0),
    stringsAsFactors     = FALSE
  )
}

#' Attach a transmission tree to a \code{diffnet} object
#'
#' Populates the \code{$transmission} slot of a \code{diffnet} with a
#' transmission tree (who-infected-whom). The resulting directed forest is the
#' canonical input to offspring-distribution analyses
#' (Lloyd-Smith \emph{et al.}, 2005) and to likelihood-based estimators of the
#' reproduction number and serial interval (White & Pagano, 2008).
#'
#' @param x A \code{diffnet} object.
#' @param tree A \code{data.frame} with at least the columns \code{date},
#'   \code{source}, \code{target}, and \code{source_exposure_date}. Columns
#'   \code{virus_id} and \code{virus} are optional. \code{source} and
#'   \code{source_exposure_date} may be \code{NA} for seed infections (roots
#'   of the tree).
#' @param pars Optional named list stored verbatim in \code{x$transmission$pars}.
#'   Useful for recording kernel parameters, seeds, etc.
#'
#' @details
#' Each row of \code{tree} represents one infection event (an edge
#' \eqn{\text{source} \to \text{target}} in the transmission tree) time-stamped
#' by \code{date}. \code{source} and \code{target} must be integer row indices
#' into \code{x} (\code{1..nnodes(x)}); \code{target} is required for every
#' row. Existing \code{$transmission} content is overwritten.
#'
#' Attaching a transmission tree promotes \code{x} to the
#' \code{\link{diffnet_epi}} subclass (\code{class(x) <- c("diffnet_epi",
#' "diffnet")}). The promotion is monotone — an already-\code{diffnet_epi}
#' input keeps its class. See \code{\link{as_diffnet_epi}} for the low-level
#' constructor.
#'
#' @return A \code{\link{diffnet_epi}} object — the input \code{x} promoted to
#'   the subclass with \code{$transmission} set to a list with components
#'   \code{tree} (a clean, ordered \code{data.frame}) and \code{pars}.
#'
#' @references
#' Lloyd-Smith, J. O., Schreiber, S. J., Kopp, P. E., & Getz, W. M. (2005).
#' Superspreading and the effect of individual variation on disease emergence.
#' \emph{Nature} 438:355-359. \doi{10.1038/nature04153}
#'
#' White, L. F., & Pagano, M. (2008). A likelihood-based method for real-time
#' estimation of the serial interval and reproductive number of an epidemic.
#' \emph{Statistics in Medicine} 27:2999-3016. \doi{10.1002/sim.3136}
#'
#' @export
#' @seealso \code{\link{new_diffnet}}
#' @author Aníbal Olivera M.
as_transmission_tree <- function(x, tree, pars = list()) {

  if (!inherits(x, "diffnet"))
    stop("-x- must be a diffnet object.")

  if (!is.data.frame(tree))
    stop("-tree- must be a data.frame.")

  required <- c("date", "source", "target", "source_exposure_date")
  missing_cols <- setdiff(required, names(tree))
  if (length(missing_cols))
    stop("-tree- is missing required column(s): ",
         paste(missing_cols, collapse = ", "), ".")

  if (anyNA(tree$target))
    stop("-tree$target- cannot contain NA values.")

  n <- nnodes(x)
  tgt <- suppressWarnings(as.integer(tree$target))
  if (anyNA(tgt) || any(tgt < 1L) || any(tgt > n))
    stop("-tree$target- must be integer indices in 1..", n, ".")

  src <- suppressWarnings(as.integer(tree$source))
  src_ok <- src[!is.na(src)]
  if (length(src_ok) && (any(src_ok < 1L) || any(src_ok > n)))
    stop("-tree$source- must be NA or an integer index in 1..", n, ".")

  out <- data.frame(
    date                 = as.integer(tree$date),
    source               = src,
    target               = tgt,
    source_exposure_date = as.integer(tree$source_exposure_date),
    stringsAsFactors     = FALSE
  )

  out$virus_id <- if (!is.null(tree$virus_id))
    as.integer(tree$virus_id) else rep(1L, nrow(out))

  out$virus <- if (!is.null(tree$virus))
    as.character(tree$virus) else rep(NA_character_, nrow(out))

  out <- out[, .transmission_cols, drop = FALSE]
  out <- out[order(out$date, out$target), , drop = FALSE]
  rownames(out) <- NULL

  # Promote to diffnet_epi (monotone — keeps class if already promoted) and
  # attach the validated tree.
  if (!inherits(x, "diffnet_epi"))
    class(x) <- c("diffnet_epi", class(x))
  x$transmission <- list(tree = out, pars = pars)
  x
}

#' Reconstruct a transmission tree from observed adoption times (M13)
#'
#' Given a dynamic contact network and per-node times of adoption, infer a
#' transmission tree by source-attribution: for each infection event, pick
#' the most plausible infector among the target's adopted neighbours at the
#' slice where the target was infected. The selection rule is the user's
#' choice (one of the bundled \code{\link{source_attribution}} kernels or a
#' user-supplied function that follows the same contract).
#'
#' This is the general-purpose primitive behind \code{rdiffnet()}'s
#' lineage-tracking (M8): the same algorithm that constructs the tree
#' during simulation also constructs it post-hoc from observed data, which
#' is what makes data products like \code{\link{epigamesDiffNet}} possible
#' without bespoke parsing code.
#'
#' @param x Either a \code{\link{diffnet}} object (graphs and \code{toa}
#'   read from its slots), or a list of adjacency matrices — one per
#'   time slice. When a list, \code{toa} must be supplied.
#' @param toa Times-of-adoption. Integer vector of length \eqn{n}
#'   (single-behaviour) or \eqn{n \times Q} integer matrix (multi-behaviour).
#'   \code{NA} marks a node that never adopted. Ignored when \code{x} is a
#'   diffnet (read from \code{x$toa} instead).
#' @param attribution The source-attribution rule. Either a string —
#'   \code{"uniform"}, \code{"weighted"}, or \code{"earliest"} (the
#'   bundled kernels) — or a function with the same signature as
#'   \code{\link{source_attribution_uniform}}.
#' @param pars Optional named list. Stored verbatim in
#'   \code{x$transmission$pars} when the result is wired into
#'   \code{\link{as_diffnet_epi}}; also forwarded to the attribution
#'   function as its \code{pars} argument.
#' @param behavior Optional character vector of length \eqn{Q} naming
#'   each diffusion process. Used to populate the \code{virus} column of
#'   the returned tree. Defaults to \code{"behavior_1"}, \code{"behavior_2"},
#'   ... when \code{NULL}.
#' @param seed Optional integer. When non-\code{NULL}, calls
#'   \code{\link{set.seed}} once before reconstruction so the stochastic
#'   attribution rules (\code{"uniform"}, \code{"weighted"}) produce a
#'   reproducible tree.
#'
#' @return A \code{data.frame} with the canonical transmission-tree
#'   schema: \code{date}, \code{source}, \code{target},
#'   \code{source_exposure_date}, \code{virus_id}, \code{virus}.
#'   One row per infection event (seeds included, with
#'   \code{source = NA}). Suitable to pass straight to
#'   \code{\link{as_transmission_tree}}.
#'
#' @details
#' For every \code{(target, virus_id)} pair with non-\code{NA} \code{toa},
#' the algorithm inspects the slice \code{x$graph[[toa[target, q]]]} and
#' picks the source among target's neighbours that were already adopted
#' (\code{toa[v, q] < toa[target, q]}). Targets with no adopted neighbour
#' at the moment of their adoption become seeds (\code{source = NA}).
#' Multi-behaviour diffnets are handled one behaviour at a time; the
#' resulting rows are concatenated.
#'
#' Under SIRS-style re-infection \code{toa[target, q]} only stores the
#' \emph{latest} infection time, so the reconstructed tree will only carry
#' one row per node per virus. To capture every entry into I, build the
#' tree at simulation time via
#' \code{rdiffnet(..., source_attribution = ...)} instead — that path
#' records each fresh adoption as it happens.
#'
#' @examples
#' set.seed(2026)
#' # Build a tiny absorbing diffnet, then reconstruct its tree post-hoc.
#' g  <- lapply(1:5, function(t) rgraph_ba(t = 4L))
#' toa <- c(1L, 2L, 3L, NA, 5L)
#' dn  <- new_diffnet(g, toa = toa, t0 = 1L, t1 = 5L)
#'
#' tree <- transmission_tree_from_events(dn, attribution = "uniform",
#'                                       seed = 2026)
#' head(tree)
#'
#' # Promote to diffnet_epi in one step:
#' dn_epi <- as_diffnet_epi(dn, attribution = "uniform")
#' is.diffnet_epi(dn_epi)
#'
#' @seealso \code{\link{source_attribution}},
#'   \code{\link{as_transmission_tree}}, \code{\link{as_diffnet_epi}}
#' @author Aníbal Olivera M.
#' @export
transmission_tree_from_events <- function(x,
                                          toa         = NULL,
                                          attribution = "uniform",
                                          pars        = list(),
                                          behavior    = NULL,
                                          seed        = NULL) {

  # 1. Coerce inputs.
  if (inherits(x, "diffnet")) {
    graphs <- x$graph
    if (is.null(toa)) toa <- x$toa
  } else if (is.list(x)) {
    graphs <- x
    if (is.null(toa))
      stop("-toa- is required when -x- is a list of graphs.")
  } else {
    stop("-x- must be a diffnet or a list of adjacency matrices.")
  }

  # 2. Normalize attribution -> function.
  attr_fn <- if (is.character(attribution)) {
    switch(attribution,
      "uniform"  = source_attribution_uniform,
      "weighted" = source_attribution_weighted,
      "earliest" = source_attribution_earliest,
      stop("Unknown -attribution-: ", attribution,
           ". Expected \"uniform\", \"weighted\", or \"earliest\".")
    )
  } else if (is.function(attribution)) {
    attribution
  } else {
    stop("-attribution- must be a function or a string.")
  }

  # 3. Normalize toa to n x Q integer matrix.
  if (is.null(dim(toa)))
    toa <- matrix(as.integer(toa), ncol = 1L)
  else
    storage.mode(toa) <- "integer"

  n <- nrow(toa)
  Q <- ncol(toa)

  # 4. Normalize behavior labels.
  if (is.null(behavior))
    behavior <- paste0("behavior_", seq_len(Q))
  else
    behavior <- as.character(behavior)
  if (length(behavior) != Q)
    stop("-behavior- must have length ", Q,
         " (one entry per column of -toa-).")

  # 5. Seed once for reproducibility of stochastic attributors.
  if (!is.null(seed)) set.seed(seed)

  # 6. Walk adopters in chronological order per behaviour.
  T_slices <- length(graphs)
  tree_rows <- list()

  for (q in seq_len(Q)) {
    adopters <- which(!is.na(toa[, q]))
    if (!length(adopters)) next
    adopters <- adopters[order(toa[adopters, q])]

    for (target in adopters) {
      t_inf <- toa[target, q]
      if (is.na(t_inf) || t_inf < 1L || t_inf > T_slices) next

      g_slice <- graphs[[t_inf]]
      row_i   <- as.vector(g_slice[target, ])
      col_i   <- as.vector(g_slice[, target])
      nbrs    <- which((row_i != 0) | (col_i != 0))
      # Keep only neighbours adopted strictly before the target.
      nbrs    <- nbrs[!is.na(toa[nbrs, q]) & toa[nbrs, q] < t_inf]

      if (length(nbrs)) {
        # Sort by ascending toa so attributors that exploit ordering
        # (e.g. source_attribution_earliest) see the same contract as
        # they do inside rdiffnet's M8 path.
        nbrs    <- nbrs[order(toa[nbrs, q])]
        weights <- pmax(row_i[nbrs], col_i[nbrs])
        src     <- attr_fn(target, nbrs, weights, t_inf, pars)
      } else {
        src <- NA_integer_
      }

      sed <- if (is.na(src)) NA_integer_ else as.integer(toa[src, q])
      tree_rows[[length(tree_rows) + 1L]] <- list(
        date                 = as.integer(t_inf),
        source               = as.integer(src),
        target               = as.integer(target),
        source_exposure_date = sed,
        virus_id             = as.integer(q),
        virus                = behavior[[q]]
      )
    }
  }

  rdiffnet_tree_rows_to_df(tree_rows)
}

#' Retrieve the transmission tree of a \code{\link{diffnet_epi}} object
#'
#' Returns the data.frame stored in \code{x$transmission$tree} for objects
#' that inherit from \code{diffnet_epi}. Plain (non-epi) diffnets do not
#' carry a tree by design; calling this function on one is an API error.
#'
#' @param x A \code{\link{diffnet_epi}} object.
#' @return A \code{data.frame} with columns \code{date}, \code{source},
#'   \code{target}, \code{source_exposure_date}, \code{virus_id}, \code{virus}.
#'   Zero rows when the epi object has been promoted but no tree attached yet.
#' @export
#' @seealso \code{\link{as_transmission_tree}}, \code{\link{as_diffnet_epi}}
#' @author Aníbal Olivera M.
transmission_tree <- function(x) {
  if (!inherits(x, "diffnet_epi"))
    stop("-x- is not a -diffnet_epi-. Use -as_transmission_tree()- or ",
         "-as_diffnet_epi()- to promote a plain diffnet first.")

  tr <- x$transmission$tree
  if (is.null(tr)) .empty_transmission_tree() else tr
}
