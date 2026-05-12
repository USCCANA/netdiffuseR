#' The \code{diffnet_epi} subclass
#'
#' \code{diffnet_epi} is an S3 subclass of \code{\link{diffnet}} carrying the
#' epidemiological extension: a \code{$transmission} slot with a
#' who-infected-whom tree and, eventually, the methods that operate on it
#' (offspring-distribution analysis, secondary attack rate, generation time,
#' survival, reproduction number, transmission-tree visualisation). A
#' \code{diffnet_epi} \emph{is a} \code{diffnet}: every method defined for the
#' base class dispatches transparently on the subclass thanks to S3
#' inheritance.
#'
#' @param x A \code{diffnet} object.
#' @param transmission Either \code{NULL} (creates an empty epi diffnet), a
#'   pre-built transmission list with the components \code{tree} and
#'   \code{pars}, or — for the data.frame entry point — pass the data.frame
#'   to \code{\link{as_transmission_tree}} instead. See examples.
#' @param pars Optional named list stored verbatim in
#'   \code{x$transmission$pars}. Only consulted when \code{transmission} is
#'   \code{NULL} or carries no \code{pars} of its own.
#' @param ... Further arguments. Accepted for compatibility with the
#'   \code{\link[base]{print}} generic; currently ignored by
#'   \code{print.diffnet_epi}.
#'
#' @return
#' \describe{
#'  \item{\code{as_diffnet_epi(x, ...)}}{A \code{diffnet_epi} object — the
#'    input \code{x} with \code{class(x) <- c("diffnet_epi", "diffnet")} and
#'    a \code{$transmission} slot.}
#'  \item{\code{is.diffnet_epi(x)}}{\code{TRUE} iff \code{x} inherits from
#'    \code{diffnet_epi}.}
#'  \item{\code{print(x)} for a \code{diffnet_epi}}{Same output as
#'    \code{\link{print.diffnet}} plus a final line summarising the
#'    transmission tree.}
#' }
#'
#' @details
#' The transmission tree is the canonical input to offspring-distribution
#' analyses (Lloyd-Smith \emph{et al.}, 2005) and likelihood-based estimators
#' of the reproduction number and serial interval (White & Pagano, 2008).
#' Attaching one to a diffnet is what turns it into an epidemiological object
#' in this package's sense — hence the dedicated subclass.
#'
#' Promotion is monotone: once a diffnet has been promoted, it stays a
#' \code{diffnet_epi}. Attaching an empty tree (\code{transmission = NULL}) is
#' allowed and useful for downstream code that wants to build the tree
#' incrementally — e.g. \code{rdiffnet()} with a future \code{source_attribution}
#' callback (M8).
#'
#' @section Class hierarchy:
#' \preformatted{
#'   class(x) -> c("diffnet_epi", "diffnet")
#' }
#' S3 dispatch tries \code{*.diffnet_epi} first, then falls back to
#' \code{*.diffnet}. \code{print.diffnet_epi} chains into
#' \code{print.diffnet} via \code{NextMethod()}.
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
#' @examples
#' set.seed(2026)
#' gr  <- lapply(1:5, function(t) rgraph_ba(t = 4L))
#' dn  <- new_diffnet(gr, toa = c(1L, 2L, NA, 3L, 5L), t0 = 1L, t1 = 5L)
#'
#' # Empty promotion (no tree yet)
#' dn_epi <- as_diffnet_epi(dn)
#' is.diffnet_epi(dn_epi)              # TRUE
#' inherits(dn_epi, "diffnet")          # also TRUE
#' nrow(transmission_tree(dn_epi))      # 0
#'
#' # Attach a tree (preferred entry point: as_transmission_tree)
#' tree <- data.frame(
#'   date = c(1L, 2L),
#'   source = c(NA_integer_, 1L),
#'   target = c(1L, 2L),
#'   source_exposure_date = c(NA_integer_, 1L)
#' )
#' dn_epi <- as_transmission_tree(dn, tree)
#' is.diffnet_epi(dn_epi)              # TRUE (promoted automatically)
#' transmission_tree(dn_epi)            # 2 rows
#'
#' @name diffnet_epi
#' @aliases diffnet_epi
#' @author Aníbal Olivera M.
NULL

#' @rdname diffnet_epi
#' @export
as_diffnet_epi <- function(x, transmission = NULL, pars = list()) {

  if (!inherits(x, "diffnet"))
    stop("-x- must be a diffnet object.")

  # Monotone promotion: prepend diffnet_epi to the class vector if absent.
  if (!inherits(x, "diffnet_epi"))
    class(x) <- c("diffnet_epi", class(x))

  if (is.null(transmission)) {
    # Empty epi diffnet (allowed by design).
    if (is.null(x$transmission))
      x$transmission <- list(tree = .empty_transmission_tree(), pars = pars)
  } else if (is.list(transmission) &&
             all(c("tree", "pars") %in% names(transmission))) {
    # Pre-built transmission list.
    x$transmission <- transmission
  } else {
    stop("-transmission- must be NULL or a list with -tree- and -pars-. ",
         "To attach a data.frame tree, use -as_transmission_tree()- instead.")
  }

  x
}

#' @rdname diffnet_epi
#' @export
is.diffnet_epi <- function(x) inherits(x, "diffnet_epi")

#' @rdname diffnet_epi
#' @export
print.diffnet_epi <- function(x, ...) {
  # Delegate the base diffnet body via NextMethod, then append one line.
  NextMethod()

  tr      <- x$transmission$tree
  n_edges <- if (is.null(tr)) 0L else nrow(tr)
  if (n_edges > 0L) {
    n_seeds <- sum(is.na(tr$source))
    n_virs  <- length(unique(tr$virus_id))
    cat(sprintf(
      "\n Transmission tree  : %d edges, %d seeds, %d virus%s",
      n_edges, n_seeds, n_virs, if (n_virs == 1L) "" else "es"
    ))
  } else {
    cat("\n Transmission tree  : empty (use -as_transmission_tree()-)")
  }

  invisible(x)
}
