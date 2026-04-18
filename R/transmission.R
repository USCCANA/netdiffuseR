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
#' @return The \code{diffnet} object \code{x} with \code{$transmission} set to
#'   a list with components \code{tree} (a clean, ordered \code{data.frame})
#'   and \code{pars}.
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

  x$transmission <- list(tree = out, pars = pars)
  x
}

#' Retrieve the transmission tree of a \code{diffnet} object
#'
#' Returns the data.frame stored in \code{x$transmission$tree}. If none has
#' been attached, an empty data.frame with the standard columns is returned.
#'
#' @param x A \code{diffnet} object.
#' @return A \code{data.frame} with columns \code{date}, \code{source},
#'   \code{target}, \code{source_exposure_date}, \code{virus_id}, \code{virus}.
#' @export
#' @seealso \code{\link{as_transmission_tree}}
transmission_tree <- function(x) {
  if (!inherits(x, "diffnet"))
    stop("-x- must be a diffnet object.")

  tr <- x$transmission$tree
  if (is.null(tr)) .empty_transmission_tree() else tr
}
