# Epidemiological metrics for diffnet / diffnet_epi (M10).
#
# Five exported metrics + a -summary.diffnet_epi- method that brings them
# together. The split between which dispatch on -diffnet- vs -diffnet_epi-:
#
#   diffnet      -> peak_prevalence, peak_time, survival_curve
#                   (need only the $status array; useful for both behaviour-
#                    diffusion and epi work)
#
#   diffnet_epi  -> secondary_attack_rate, generation_time
#                   (need the $transmission tree)
#
# The user can invoke any of the first three on a plain diffnet; the latter
# two error on a plain diffnet pointing to as_transmission_tree() /
# as_diffnet_epi(). summary.diffnet_epi prints the base diffnet block plus an
# epi block that calls all five.

# ----------------------------------------------------------------------------
# peak_prevalence / peak_time
# ----------------------------------------------------------------------------

#' Peak prevalence and time of peak in a diffnet
#'
#' \code{peak_prevalence(x)} returns the highest fraction of adopted nodes
#' observed in \code{x$status} across all time slices. \code{peak_time(x)}
#' returns the time period at which the peak is reached.
#'
#' For multi-behaviour diffnets both functions return a named numeric vector
#' of length \eqn{Q}, one entry per behaviour.
#'
#' @param x A \code{\link{diffnet}} (or any subclass, such as
#'   \code{\link{diffnet_epi}}).
#' @param ... Currently ignored.
#'
#' @return Numeric scalar (single behaviour) or named numeric vector
#'   (multi-behaviour). \code{peak_prevalence} is in \eqn{[0, 1]};
#'   \code{peak_time} is the time-period label (as integer).
#'
#' @examples
#' set.seed(2026)
#' dn <- rdiffnet(n = 50, t = 8, seed.graph = "small-world",
#'                seed.p.adopt = 0.05, stop.no.diff = FALSE)
#' peak_prevalence(dn)
#' peak_time(dn)
#'
#' @name peak_prevalence
#' @author Aníbal Olivera M.
NULL

#' @rdname peak_prevalence
#' @export
peak_prevalence <- function(x, ...) UseMethod("peak_prevalence")

#' @rdname peak_prevalence
#' @export
peak_prevalence.diffnet <- function(x, ...) {
  st <- x$status
  if (is.list(st)) {
    res <- vapply(st, function(s) max(colSums(s)) / nrow(s), numeric(1L))
    names(res) <- vapply(x$meta$behavior, as.character, character(1L))
    return(res)
  }
  max(colSums(st)) / nrow(st)
}

#' @rdname peak_prevalence
#' @export
peak_time <- function(x, ...) UseMethod("peak_time")

#' @rdname peak_prevalence
#' @export
peak_time.diffnet <- function(x, ...) {
  st   <- x$status
  pers <- x$meta$pers
  if (is.list(st)) {
    res <- vapply(st, function(s) as.integer(pers[which.max(colSums(s))]),
                  integer(1L))
    names(res) <- vapply(x$meta$behavior, as.character, character(1L))
    return(res)
  }
  as.integer(pers[which.max(colSums(st))])
}

# ----------------------------------------------------------------------------
# survival_curve
# ----------------------------------------------------------------------------

#' Kaplan-Meier-style survival curve for a diffnet
#'
#' For each adopted node, \code{survival_curve(x)} computes the duration in
#' the adopted state (\code{tod(x) - toa(x)}; right-censored at \code{T} for
#' nodes that never recover) and assembles a Kaplan-Meier-style survival
#' table.
#'
#' @param x A \code{\link{diffnet}} object.
#' @param ... Currently ignored.
#'
#' @return A \code{data.frame} (with extra class \code{netdiffuseR_survival})
#'   carrying columns \code{time}, \code{n_at_risk}, \code{n_recovered}, and
#'   \code{survival}. For multi-behaviour diffnets an additional
#'   \code{virus_id} column tags the behaviour. Printing summarises the
#'   curve; standard data.frame subscripting works on the underlying rows.
#'
#' @examples
#' set.seed(2026)
#' dn <- rdiffnet(n = 50, t = 8, seed.graph = "small-world",
#'                seed.p.adopt = 0.10, stop.no.diff = FALSE,
#'                disadopt = disadoptmech_random(prob = 0.15))
#' s <- survival_curve(dn)
#' s                           # prints summary
#' as.data.frame(s)             # full data.frame
#'
#' @name survival_curve
#' @author Aníbal Olivera M.
NULL

#' @rdname survival_curve
#' @export
survival_curve <- function(x, ...) UseMethod("survival_curve")

#' @rdname survival_curve
#' @export
survival_curve.diffnet <- function(x, ...) {
  toa_v <- x$toa
  tod_v <- tod(x)
  T_end <- x$meta$nper

  if (is.null(dim(toa_v))) {
    out <- .survival_one_behaviour(toa_v, tod_v, T_end)
  } else {
    Q <- ncol(toa_v)
    parts <- lapply(seq_len(Q), function(q) {
      df <- .survival_one_behaviour(toa_v[, q], tod_v[, q], T_end)
      df$virus_id <- q
      df[, c("virus_id", "time", "n_at_risk", "n_recovered", "survival")]
    })
    out <- do.call(rbind, parts)
  }

  rownames(out) <- NULL
  structure(out, class = c("netdiffuseR_survival", "data.frame"))
}

# Internal: KM table for one behaviour.
.survival_one_behaviour <- function(toa, tod, T_end) {
  adopters <- which(!is.na(toa))
  if (!length(adopters)) {
    return(data.frame(time = integer(0), n_at_risk = integer(0),
                      n_recovered = integer(0), survival = numeric(0),
                      stringsAsFactors = FALSE))
  }

  # duration in adopted state, censored at T_end - toa + 1 when absorbing
  duration  <- ifelse(is.na(tod[adopters]),
                      T_end - toa[adopters] + 1L,
                      tod[adopters] - toa[adopters])
  recovered <- !is.na(tod[adopters])

  unique_durs <- sort(unique(duration))
  surv      <- 1
  remaining <- length(adopters)

  rows <- vector("list", length(unique_durs))
  for (k in seq_along(unique_durs)) {
    d   <- unique_durs[k]
    rec <- sum(duration == d & recovered)
    cen <- sum(duration == d & !recovered)
    if (rec > 0 && remaining > 0) surv <- surv * (1 - rec / remaining)
    rows[[k]] <- list(time = as.integer(d),
                      n_at_risk = as.integer(remaining),
                      n_recovered = as.integer(rec),
                      survival = as.numeric(surv))
    remaining <- remaining - rec - cen
  }

  data.frame(
    time        = vapply(rows, `[[`, integer(1L), "time"),
    n_at_risk   = vapply(rows, `[[`, integer(1L), "n_at_risk"),
    n_recovered = vapply(rows, `[[`, integer(1L), "n_recovered"),
    survival    = vapply(rows, `[[`, numeric(1L), "survival"),
    stringsAsFactors = FALSE
  )
}

#' @rdname survival_curve
#' @export
print.netdiffuseR_survival <- function(x, ...) {
  cat("Survival curve (Kaplan-Meier-style)\n")
  if (!nrow(x)) {
    cat(" Empty -- no adopters in this diffnet.\n")
    return(invisible(x))
  }
  has_virus <- "virus_id" %in% names(x)
  vids <- if (has_virus) unique(x$virus_id) else 1L
  for (v in vids) {
    sub <- if (has_virus) x[x$virus_id == v, , drop = FALSE] else x
    n_events <- sum(sub$n_recovered)
    surv_at_end <- sub$survival[nrow(sub)]
    # Median: first time where survival <= 0.5
    median_t <- if (any(sub$survival <= 0.5))
      sub$time[which(sub$survival <= 0.5)[1L]] else NA_integer_
    tag <- if (has_virus) sprintf(" [behaviour %d]", v) else ""
    cat(sprintf(" N events%s : %d recoveries across %d distinct durations\n",
                tag, n_events, nrow(sub)))
    cat(sprintf(" Median survival : %s\n",
                if (is.na(median_t)) "not reached"
                else as.character(median_t)))
    cat(sprintf(" Final survival (end of horizon) : %.3f\n", surv_at_end))
  }
  cat(" -> use as.data.frame(.) for the per-time table.\n")
  invisible(x)
}

# ----------------------------------------------------------------------------
# secondary_attack_rate (diffnet_epi only)
# ----------------------------------------------------------------------------

#' Secondary attack rate from a transmission tree
#'
#' For each infection event recorded in \code{$transmission$tree},
#' \code{secondary_attack_rate(x)} reports the number of secondary infections
#' caused by that event and the number of contacts the infector had in the
#' contact network at the slice corresponding to \code{source_exposure_date}
#' (the infector's own infection date). The per-event rate is
#' \code{n_secondary / n_contacts}; the aggregate (printed by default) is
#' \code{sum(n_secondary) / sum(n_contacts)}.
#'
#' Under absorbing diffusion each \code{(source, virus_id)} has exactly one
#' \code{source_exposure_date}, so the per-event keying collapses to the
#' classic per-source rollup. Under SIRS-style re-infection (a node enters
#' state I multiple times for the same virus), each infection-life of the
#' source is its own row, matching the convention used by epiworldR for
#' tree-derived metrics.
#'
#' Under SIRS the same \code{(source, target)} pair can transmit multiple
#' times during the source's infection life (the target disadopts and gets
#' re-infected by the same source). Each such transmission is a distinct
#' row in the tree and contributes to \code{n_secondary} for that
#' source-event, while \code{n_contacts} is fixed at the source's
#' neighbourhood size at \code{source_exposure_date}. Consequently the
#' per-event \code{sar} may exceed 1 (it is no longer a probability of
#' transmission but a count of transmissions per contact). The aggregate
#' \code{attr(sar, "global")} retains its sum-over-sum interpretation.
#'
#' @param x A \code{\link{diffnet_epi}} object.
#' @param ... Currently ignored.
#'
#' @return A \code{data.frame} (with extra class \code{netdiffuseR_sar})
#'   carrying columns \code{source}, \code{virus_id},
#'   \code{source_exposure_date}, \code{n_secondary}, \code{n_contacts},
#'   \code{sar}. Printing shows the aggregate scalar; the per-event rows
#'   are exposed via standard data.frame subscripting. The aggregate is
#'   also stored as \code{attr(., "global")}.
#'
#' @examples
#' set.seed(2026)
#' dn <- rdiffnet(n = 40, t = 8, seed.graph = "small-world",
#'                seed.p.adopt = 0.10, stop.no.diff = FALSE,
#'                source_attribution = source_attribution_uniform)
#' sar <- secondary_attack_rate(dn)
#' sar                                  # aggregate print
#' as.data.frame(sar)                    # per-source breakdown
#' attr(sar, "global")                   # aggregate scalar
#'
#' @name secondary_attack_rate
#' @author Aníbal Olivera M.
NULL

#' @rdname secondary_attack_rate
#' @export
secondary_attack_rate <- function(x, ...) UseMethod("secondary_attack_rate")

#' @rdname secondary_attack_rate
#' @export
secondary_attack_rate.default <- function(x, ...) {
  stop("-secondary_attack_rate()- requires a -diffnet_epi-. ",
       "Use -as_transmission_tree()- or -as_diffnet_epi()- first.")
}

#' @rdname secondary_attack_rate
#' @export
secondary_attack_rate.diffnet_epi <- function(x, ...) {

  tr <- transmission_tree(x)
  tr_edges <- tr[!is.na(tr$source), , drop = FALSE]

  if (!nrow(tr_edges)) {
    out <- data.frame(
      source = integer(0), virus_id = integer(0),
      source_exposure_date = integer(0),
      n_secondary = integer(0), n_contacts = integer(0),
      sar = numeric(0), stringsAsFactors = FALSE
    )
    return(structure(out, global = NA_real_,
                     class = c("netdiffuseR_sar", "data.frame")))
  }

  # M12.2: aggregate secondaries per *infection event* of the source,
  # keyed by (source, virus_id, source_exposure_date). Under absorbing
  # diffusion each source has one exposure_date and this collapses to
  # the M12 per-source rollup; under SIRS-style re-infection each
  # infection-life of the source is its own row.
  key <- paste(tr_edges$source, tr_edges$virus_id,
               tr_edges$source_exposure_date, sep = "::")
  agg <- tapply(seq_len(nrow(tr_edges)), key, length)
  parts <- strsplit(names(agg), "::", fixed = TRUE)
  src   <- as.integer(vapply(parts, `[`, character(1L), 1L))
  vid   <- as.integer(vapply(parts, `[`, character(1L), 2L))
  sed   <- as.integer(vapply(parts, `[`, character(1L), 3L))
  n_sec <- as.integer(agg)

  # Contacts at the slice corresponding to *this* infection event of the
  # source, i.e. graph[[source_exposure_date]] (not graph[[toa[source]]],
  # which under re-infection only carries the latest exposure date).
  T_slices <- length(x$graph)
  n_con <- vapply(seq_along(src), function(i) {
    t_inf <- sed[i]
    if (is.na(t_inf) || t_inf < 1L || t_inf > T_slices) return(0L)
    g <- x$graph[[t_inf]]
    s <- src[i]
    as.integer(sum((g[s, ] != 0) | (g[, s] != 0)))
  }, integer(1L))

  per_event <- data.frame(
    source               = src,
    virus_id             = vid,
    source_exposure_date = sed,
    n_secondary          = n_sec,
    n_contacts           = n_con,
    sar                  = ifelse(n_con > 0, n_sec / n_con, NA_real_),
    stringsAsFactors = FALSE
  )
  per_event <- per_event[order(per_event$virus_id, per_event$source,
                                per_event$source_exposure_date), ,
                          drop = FALSE]
  rownames(per_event) <- NULL

  total_s <- sum(per_event$n_secondary)
  total_c <- sum(per_event$n_contacts)
  global  <- if (total_c > 0) total_s / total_c else NA_real_

  structure(per_event, global = global,
            class = c("netdiffuseR_sar", "data.frame"))
}

#' @rdname secondary_attack_rate
#' @export
print.netdiffuseR_sar <- function(x, ...) {
  cat("Secondary Attack Rate\n")
  cat(sprintf(" Aggregate (sum of secondaries / sum of contacts) : %.3f\n",
              attr(x, "global")))
  n_events  <- nrow(x)
  n_sources <- length(unique(paste(x$source, x$virus_id, sep = "::")))
  if (n_events == n_sources) {
    cat(sprintf(" Based on %d infector%s in the transmission tree.\n",
                n_events, if (n_events == 1L) "" else "s"))
  } else {
    cat(sprintf(" Based on %d infection event%s from %d distinct infector%s\n",
                n_events, if (n_events == 1L) "" else "s",
                n_sources, if (n_sources == 1L) "" else "s"))
    cat("  (some infectors re-entered I and transmitted in more than one life).\n")
  }
  cat(" -> use as.data.frame(.) or standard subscripting for the\n")
  cat("    per-event breakdown.\n")
  invisible(x)
}

# ----------------------------------------------------------------------------
# generation_time (diffnet_epi only)
# ----------------------------------------------------------------------------

#' Generation time per edge of a transmission tree
#'
#' For each edge \eqn{(source, target)} in \code{$transmission$tree},
#' \code{generation_time(x)} computes \code{date - source_exposure_date},
#' the time between the infector's adoption and its infectee's. Seed
#' rows (\code{source == NA}) are dropped.
#'
#' @param x A \code{\link{diffnet_epi}} object.
#' @param ... Currently ignored.
#'
#' @return A \code{data.frame} (with extra class
#'   \code{netdiffuseR_generation_time}) carrying the original tree columns
#'   plus \code{gen_time}. Printing shows a distributional summary
#'   (\emph{N}, mean, sd, median, IQR, range); the per-edge rows are exposed
#'   via standard data.frame subscripting.
#'
#' @examples
#' set.seed(2026)
#' dn <- rdiffnet(n = 40, t = 8, seed.graph = "small-world",
#'                seed.p.adopt = 0.10, stop.no.diff = FALSE,
#'                source_attribution = source_attribution_uniform)
#' gt <- generation_time(dn)
#' gt                                   # summary print
#' as.data.frame(gt)                    # per-edge rows
#' mean(gt$gen_time)
#'
#' @name generation_time
#' @author Aníbal Olivera M.
NULL

#' @rdname generation_time
#' @export
generation_time <- function(x, ...) UseMethod("generation_time")

#' @rdname generation_time
#' @export
generation_time.default <- function(x, ...) {
  stop("-generation_time()- requires a -diffnet_epi-. ",
       "Use -as_transmission_tree()- or -as_diffnet_epi()- first.")
}

#' @rdname generation_time
#' @export
generation_time.diffnet_epi <- function(x, ...) {
  tr <- transmission_tree(x)
  tr <- tr[!is.na(tr$source) & !is.na(tr$source_exposure_date), , drop = FALSE]
  tr$gen_time <- as.integer(tr$date - tr$source_exposure_date)
  rownames(tr) <- NULL
  structure(tr, class = c("netdiffuseR_generation_time", "data.frame"))
}

#' @rdname generation_time
#' @export
print.netdiffuseR_generation_time <- function(x, ...) {
  cat("Generation time distribution\n")
  if (!nrow(x)) {
    cat(" Empty -- no non-seed edges in the transmission tree.\n")
    return(invisible(x))
  }
  g <- x$gen_time
  q <- stats::quantile(g, c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
  cat(sprintf(" N edges : %d\n", length(g)))
  cat(sprintf(" Mean    : %.2f  (sd %.2f)\n", mean(g), stats::sd(g)))
  cat(sprintf(" Median  : %.1f  (IQR %.1f - %.1f)\n", q[2], q[1], q[3]))
  cat(sprintf(" Range   : %d - %d\n", min(g), max(g)))
  cat(" -> use as.data.frame(.) for the per-edge table.\n")
  invisible(x)
}

# ----------------------------------------------------------------------------
# repr_number (diffnet_epi only)
# ----------------------------------------------------------------------------

#' Empirical reproduction number from a transmission tree
#'
#' For every infection event in \code{$transmission$tree},
#' \code{repr_number(x)} counts the number of secondary cases it caused
#' (its offspring count, \eqn{\nu_i} in Lloyd-Smith \emph{et al.}, 2005)
#' and reports the mean across cases as the empirical reproduction
#' number. Cases that did not transmit further (terminal cases) count
#' as zero in the denominator; seeds are included.
#'
#' A case is one entry into state I, keyed by
#' \code{(node, virus_id, exposure_date)}. Under absorbing diffusion
#' (the classic netdiffuseR regime) each \code{(node, virus_id)} has
#' exactly one \code{exposure_date}, so the 3-D key collapses to the
#' familiar per-node rollup. Under SIRS-style re-infection (a node
#' enters \eqn{I} multiple times for the same virus), each
#' infection-life is its own case with its own offspring tally. This
#' matches the convention used by epiworldR's
#' \code{get_reproductive_number()} and the Lloyd-Smith framework.
#'
#' @param x A \code{\link{diffnet_epi}} object.
#' @param ... Currently ignored.
#'
#' @return A \code{data.frame} (with extra class \code{netdiffuseR_repr})
#'   carrying columns \code{node}, \code{virus_id},
#'   \code{exposure_date}, \code{n_offspring}. Printing shows the
#'   aggregate reproduction number (mean offspring), plus SD and range;
#'   the per-case rows are exposed via standard data.frame subscripting.
#'   The aggregate is also stored as \code{attr(., "global")}. A
#'   \code{plot} method renders the offspring distribution as a
#'   barplot.
#'
#' @details
#' The empirical reproduction number is defined as the mean offspring
#' count across all observed cases:
#'
#' \deqn{%
#' R = \frac{1}{N}\sum_{i \in \mathrm{cases}} \nu_i %
#' }{%
#' R = (1/N) * sum_i nu_i %
#' }
#'
#' where \eqn{N} is the total number of infected cases (seeds + secondary)
#' in the tree and \eqn{\nu_i} is the number of times case \eqn{i} appears
#' as a \code{source} in the tree. Terminal cases (\eqn{\nu_i = 0}) are
#' included in the denominator, so \eqn{R} is the true mean offspring,
#' not the mean among transmitters only.
#'
#' For trees built from observational data (Epigames / contact tracing),
#' \eqn{R} matches the standard tree-based reproduction-number estimator.
#' For trees produced by \code{rdiffnet()} with \code{source_attribution},
#' the value depends on the attribution policy: \code{_uniform},
#' \code{_weighted}, and \code{_earliest} will produce different empirical
#' \eqn{R} on the same simulation, since they distribute observed
#' adoptions across different infectors.
#'
#' @references
#' Lloyd-Smith, J. O., Schreiber, S. J., Kopp, P. E., & Getz, W. M. (2005).
#' Superspreading and the effect of individual variation on disease emergence.
#' \emph{Nature} 438:355-359. \doi{10.1038/nature04153}
#'
#' @examples
#' set.seed(2026)
#' dn <- rdiffnet(n = 40, t = 6, seed.graph = "small-world",
#'                seed.p.adopt = 0.10, stop.no.diff = FALSE,
#'                source_attribution = source_attribution_uniform)
#' R <- repr_number(dn)
#' R                                 # aggregate print: mean / SD / range
#' as.data.frame(R)                  # per-case offspring counts
#' attr(R, "global")                 # the scalar R
#' \dontrun{
#' plot(R)                            # offspring distribution barplot
#' }
#'
#' # SIRS-style: a disadopt function lets nodes re-enter I, and every
#' # re-infection is recorded as its own case in the returned frame.
#' \dontrun{
#' disadopt_30 <- function(expo, cumadopt, time) {
#'   q_max <- dim(cumadopt)[3]; res <- vector("list", q_max)
#'   for (q in seq_len(q_max)) {
#'     adopters <- which(cumadopt[, time, q] == 1L)
#'     res[[q]] <- if (length(adopters))
#'       sample(adopters, ceiling(0.30 * length(adopters))) else integer()
#'   }
#'   res
#' }
#' set.seed(2026)
#' dn_sirs <- rdiffnet(n = 60, t = 10, seed.graph = "small-world",
#'                     seed.p.adopt = 0.15, stop.no.diff = FALSE,
#'                     disadopt = disadopt_30,
#'                     source_attribution = source_attribution_uniform)
#' R_sirs <- repr_number(dn_sirs)
#' table(table(paste(R_sirs$node, R_sirs$virus_id))) # nodes by # of lives
#' }
#'
#' @name repr_number
#' @author Aníbal Olivera M.
NULL

#' @rdname repr_number
#' @export
repr_number <- function(x, ...) UseMethod("repr_number")

#' @rdname repr_number
#' @export
repr_number.default <- function(x, ...) {
  stop("-repr_number()- requires a -diffnet_epi-. ",
       "Use -as_transmission_tree()- or -as_diffnet_epi()- first.")
}

#' @rdname repr_number
#' @export
repr_number.diffnet_epi <- function(x, ...) {

  tr <- transmission_tree(x)

  if (!nrow(tr)) {
    out <- data.frame(
      node          = integer(0),
      virus_id      = integer(0),
      exposure_date = integer(0),
      n_offspring   = integer(0),
      stringsAsFactors = FALSE
    )
    return(structure(out, global = NA_real_,
                     class = c("netdiffuseR_repr", "data.frame")))
  }

  # M12.2: cases are keyed per infection event = unique
  # (target, virus_id, date). Under single-adoption each (target, virus_id)
  # has exactly one date, so this collapses to the M12 2-D keying for
  # absorbing diffusions; under SIRS-style re-infection each entry-to-I
  # is its own case (Lloyd-Smith / epiworldR convention).
  cases <- unique(tr[, c("target", "virus_id", "date"), drop = FALSE])
  names(cases) <- c("node", "virus_id", "exposure_date")
  key_cases <- paste(cases$node, cases$virus_id, cases$exposure_date,
                     sep = "::")

  # Count offspring per source-event: a row in the tree contributes +1 to
  # the case that infected someone *at that source_exposure_date*.
  src_rows <- tr[!is.na(tr$source),
                 c("source", "virus_id", "source_exposure_date"),
                 drop = FALSE]
  if (nrow(src_rows)) {
    key_src     <- paste(src_rows$source, src_rows$virus_id,
                         src_rows$source_exposure_date, sep = "::")
    src_tab     <- table(factor(key_src, levels = key_cases))
    n_offspring <- as.integer(src_tab)
  } else {
    n_offspring <- rep(0L, nrow(cases))
  }

  cases$n_offspring <- n_offspring
  cases             <- cases[order(cases$virus_id, cases$node,
                                   cases$exposure_date), ,
                              drop = FALSE]
  rownames(cases)   <- NULL

  global <- mean(cases$n_offspring)

  structure(cases, global = global,
            class = c("netdiffuseR_repr", "data.frame"))
}

#' @rdname repr_number
#' @export
print.netdiffuseR_repr <- function(x, ...) {
  cat("Reproduction number (empirical, from transmission tree)\n")
  if (!nrow(x)) {
    cat(" Empty -- no cases in the transmission tree.\n")
    return(invisible(x))
  }
  nv <- length(unique(x$virus_id))
  if (nv > 1L) {
    cat(sprintf(" Aggregate over %d diffusions (pooled).\n", nv))
    cat(sprintf(" Mean offspring (R) : %.3f\n", attr(x, "global")))
    if (nrow(x) > 1L)
      cat(sprintf(" SD                 : %.3f\n", stats::sd(x$n_offspring)))
    cat(sprintf(" Range              : %d - %d\n",
                min(x$n_offspring), max(x$n_offspring)))
    cat(sprintf(" Based on %d cases across %d diffusions.\n", nrow(x), nv))
    cat(" Per-diffusion R:\n")
    for (v in sort(unique(x$virus_id))) {
      sub <- x$n_offspring[x$virus_id == v]
      cat(sprintf("  diffusion %s: R = %.3f  (n = %d)\n",
                  format(v), mean(sub), length(sub)))
    }
    cat(" -> use as.data.frame(.) for the per-case offspring count,\n")
    cat("    subset by $virus_id for per-diffusion rows,\n")
    cat("    or plot(.) for the offspring distribution.\n")
  } else {
    cat(sprintf(" Mean offspring (R) : %.3f\n", attr(x, "global")))
    if (nrow(x) > 1L)
      cat(sprintf(" SD                 : %.3f\n", stats::sd(x$n_offspring)))
    cat(sprintf(" Range              : %d - %d\n",
                min(x$n_offspring), max(x$n_offspring)))
    cat(sprintf(" Based on %d case%s in the transmission tree.\n",
                nrow(x), if (nrow(x) == 1L) "" else "s"))
    cat(" -> use as.data.frame(.) for the per-case offspring count,\n")
    cat("    or plot(.) for the offspring distribution.\n")
  }
  invisible(x)
}

#' @rdname repr_number
#' @param y Unused. Present for S3 consistency with \code{\link[graphics]{plot}}.
#' @param main Plot title. When \code{NULL} (default), a sensible title is
#'   chosen automatically: it includes "pooled over k diffusions" when the
#'   tree carries multiple diffusion processes (i.e., multiple
#'   \code{virus_id} values), otherwise just "Offspring distribution".
#' @param xlab,ylab Axis labels forwarded to \code{\link[graphics]{barplot}}.
#' @export
plot.netdiffuseR_repr <- function(x, y = NULL,
                                  main = NULL,
                                  xlab = "Number of offspring (secondary cases)",
                                  ylab = "Number of cases",
                                  ...) {
  if (!nrow(x)) {
    if (is.null(main)) main <- "Offspring distribution"
    graphics::plot.new()
    graphics::title(main = main, sub = "Empty transmission tree")
    return(invisible(x))
  }
  nv <- length(unique(x$virus_id))
  if (is.null(main)) {
    main <- if (nv > 1L)
      sprintf("Offspring distribution (pooled over %d diffusions)", nv)
    else
      "Offspring distribution"
  }
  sub <- if (nv > 1L)
    sprintf("R = %.3f across %d cases / %d diffusions",
            attr(x, "global"), nrow(x), nv)
  else
    sprintf("R = %.3f across %d cases", attr(x, "global"), nrow(x))

  k   <- max(x$n_offspring)
  tab <- table(factor(x$n_offspring, levels = 0:k))
  graphics::barplot(as.numeric(tab),
                    names.arg = names(tab),
                    main = main, sub = sub,
                    xlab = xlab, ylab = ylab, ...)
  invisible(x)
}

# ----------------------------------------------------------------------------
# summary.diffnet_epi
# ----------------------------------------------------------------------------

#' Summary method for \code{diffnet_epi} objects
#'
#' Extends \code{\link[=summary.diffnet]{summary.diffnet}} with an
#' epidemiological block: peak prevalence + peak time, secondary attack rate
#' (aggregate), and generation time (summary stats). The base diffnet block
#' is printed first via \code{NextMethod()}; the epi block follows.
#'
#' @param object A \code{\link{diffnet_epi}} object.
#' @param ... Forwarded to \code{summary.diffnet}.
#'
#' @return Invisibly, the object returned by \code{summary.diffnet}. The
#'   epi block is printed as a side effect.
#' @export
#' @author Aníbal Olivera M.
summary.diffnet_epi <- function(object, ...) {
  base_summary <- NextMethod()

  cat("\n Epidemiological metrics ----------------------------------\n")

  # Peak prevalence and peak time (work via the base diffnet method)
  pp <- peak_prevalence(object)
  pt <- peak_time(object)
  if (length(pp) == 1L) {
    cat(sprintf(" Peak prevalence : %.3f at t = %d\n", pp, pt))
  } else {
    for (q in seq_along(pp)) {
      cat(sprintf(" Peak prevalence [%s] : %.3f at t = %d\n",
                  names(pp)[q], pp[q], pt[q]))
    }
  }

  # SAR (diffnet_epi-only; safe to call here)
  sar <- secondary_attack_rate(object)
  cat(sprintf(" Secondary Attack Rate (aggregate) : %.3f (n infectors = %d)\n",
              attr(sar, "global"), nrow(sar)))

  # Generation time (diffnet_epi-only)
  gt <- generation_time(object)
  if (nrow(gt) > 0L) {
    cat(sprintf(
      " Generation time : mean %.2f, median %.1f (N edges = %d)\n",
      mean(gt$gen_time), stats::median(gt$gen_time), nrow(gt)
    ))
  } else {
    cat(" Generation time : empty (seed-only tree)\n")
  }

  # Survival curve (works on plain diffnet too; meaningful only when
  # there is at least one recovery event in $status)
  sc <- survival_curve(object)
  if (nrow(sc) > 0L) {
    n_rec   <- sum(sc$n_recovered)
    surv_end <- sc$survival[nrow(sc)]
    if (n_rec > 0L) {
      med_ix <- which(sc$survival <= 0.5)[1L]
      med_t  <- if (length(med_ix) && !is.na(med_ix)) sc$time[med_ix] else NA_integer_
      cat(sprintf(
        " Survival curve  : %d recoveries; final survival %.3f%s\n",
        n_rec, surv_end,
        if (is.na(med_t)) "" else sprintf("; median t = %s", as.character(med_t))
      ))
    } else {
      cat(" Survival curve  : flat at 1 (no recoveries observed)\n")
    }
  } else {
    cat(" Survival curve  : empty (no adopters)\n")
  }

  invisible(base_summary)
}
