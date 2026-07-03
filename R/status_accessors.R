#' Accessors for adoption / disadoption times in a \code{diffnet}
#'
#' Mirror accessors over the canonical \code{$status} slot. \code{toa} and
#' \code{tod} return the \emph{first} adoption / first recovery time per node
#' (per behaviour for multi-behaviour diffnets) — same shape as the legacy
#' \code{$toa} field. \code{toa_all} and \code{tod_all} return long-format
#' \code{data.frame}s capturing every event in the multi-cycle history.
#'
#' @param x A \code{diffnet} object.
#'
#' @return
#' \describe{
#'  \item{\code{toa(x)}}{Integer vector of length \eqn{n} (single-behaviour)
#'    or \eqn{n \times Q} integer matrix (multi-behaviour). \code{NA} when
#'    the node never adopted. Equivalent to \code{x$toa}.}
#'  \item{\code{tod(x)}}{Same shape as \code{toa(x)}. First time after
#'    \code{toa[i, q]} when \code{$status[i, t, q]} flips back to 0.
#'    \code{NA} when the node never recovered (absorbing).}
#'  \item{\code{toa_all(x)}}{\code{data.frame} with columns \code{node},
#'    \code{behavior}, \code{episode}, \code{time}. One row per fresh
#'    adoption event in \code{$status}.}
#'  \item{\code{tod_all(x)}}{\code{data.frame} with columns \code{node},
#'    \code{behavior}, \code{episode}, \code{time}. One row per recovery
#'    event in \code{$status}.}
#' }
#'
#' @details
#' For an absorbing single-cycle diffnet \code{tod(x)} is \code{NA} for every
#' node (no recoveries) and \code{tod_all(x)} returns a zero-row data.frame.
#' For a multi-cycle diffnet, \code{tod(x)} reports only the first recovery
#' per node-behaviour as a summary; use \code{tod_all(x)} for the full
#' history.
#'
#' @examples
#' set.seed(2026)
#' g  <- rgraph_er(n = 10, t = 1, p = 0.4)
#' dn <- rdiffnet(seed.graph = g, t = 6, seed.p.adopt = 0.2,
#'                stop.no.diff = FALSE)
#'
#' toa(dn)        # first adoption time per node (same as dn$toa)
#' tod(dn)        # first recovery — all NA for an absorbing diffnet
#' toa_all(dn)    # one row per fresh adoption event
#' tod_all(dn)    # zero-row data.frame for an absorbing diffnet
#'
#' @author Aníbal Olivera M.
#' @name status_accessors
NULL

# ----------------------------------------------------------------------------
# toa / tod : simple summaries (same shape as the $toa slot)
# ----------------------------------------------------------------------------

#' @rdname status_accessors
#' @export
toa <- function(x) UseMethod("toa")

#' @rdname status_accessors
#' @export
toa.diffnet <- function(x) x$toa

#' @rdname status_accessors
#' @export
tod <- function(x) UseMethod("tod")

#' @rdname status_accessors
#' @export
tod.diffnet <- function(x) {
  s <- x$status
  if (is.null(s)) {
    # Defensive — every diffnet built post-status-refactor has a $status
    # slot, but keep this safe against any older object that might still
    # be in scope.
    return(rep(NA_integer_, length(x$toa)))
  }

  if (is.list(s)) {
    Q <- length(s)
    n <- nrow(s[[1L]])
    out <- matrix(NA_integer_, n, Q)
    for (q in seq_len(Q)) out[, q] <- .first_recovery(s[[q]], x$toa[, q])
    rownames(out) <- rownames(x$toa)
    return(out)
  }

  res <- .first_recovery(s, x$toa)
  if (length(names(x$toa))) names(res) <- names(x$toa)
  res
}

# Internal: first time t > toa[i] when status[i, t] == 0; NA otherwise.
.first_recovery <- function(status_q, toa_q) {
  n <- nrow(status_q)
  T <- ncol(status_q)
  vapply(seq_len(n), function(i) {
    ti <- toa_q[i]
    if (is.na(ti) || ti >= T) return(NA_integer_)
    after <- which(status_q[i, (ti + 1L):T] == 0L)
    if (length(after)) ti + as.integer(after[1L]) else NA_integer_
  }, integer(1L))
}

# ----------------------------------------------------------------------------
# toa_all / tod_all : long-format multi-cycle history
# ----------------------------------------------------------------------------

#' @rdname status_accessors
#' @export
toa_all <- function(x) UseMethod("toa_all")

#' @rdname status_accessors
#' @export
toa_all.diffnet <- function(x) .events_long(x, kind = "adopt")

#' @rdname status_accessors
#' @export
tod_all <- function(x) UseMethod("tod_all")

#' @rdname status_accessors
#' @export
tod_all.diffnet <- function(x) .events_long(x, kind = "recover")

# Internal: walk $status (single matrix or list of matrices) and emit a
# long-format data.frame with one row per "fresh adoption" (kind = "adopt")
# or per "fresh recovery" (kind = "recover") event. The episode column is
# 1-indexed within each (node, behavior) pair, in time order.
.events_long <- function(x, kind = c("adopt", "recover")) {
  kind <- match.arg(kind)
  s <- x$status

  if (is.null(s)) {
    return(data.frame(
      node     = integer(0),
      behavior = integer(0),
      episode  = integer(0),
      time     = integer(0)
    ))
  }

  status_list <- if (is.list(s)) s else list(s)

  rows <- list()
  for (q in seq_along(status_list)) {
    sq <- status_list[[q]]
    n  <- nrow(sq); T <- ncol(sq)
    times_label <- as.integer(colnames(sq))
    if (length(times_label) != T) times_label <- seq_len(T)

    for (i in seq_len(n)) {
      v <- as.integer(sq[i, ])
      # diffs[t] = v[t] - v[t - 1]; +1 = fresh adoption, -1 = fresh recovery.
      diffs <- c(v[1L], diff(v))
      events_t <- if (kind == "adopt") which(diffs == 1L) else which(diffs == -1L)
      if (length(events_t)) {
        rows[[length(rows) + 1L]] <- data.frame(
          node     = rep.int(i, length(events_t)),
          behavior = rep.int(q, length(events_t)),
          episode  = seq_along(events_t),
          time     = times_label[events_t],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(rows)) {
    return(data.frame(
      node     = integer(0),
      behavior = integer(0),
      episode  = integer(0),
      time     = integer(0)
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$node, out$behavior, out$episode), , drop = FALSE]
}
