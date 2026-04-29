#' Disadoption mechanisms for \code{rdiffnet}
#'
#' A family of factories that build per-step disadoption rules. Each
#' factory takes its parameters and returns a closure that satisfies the
#' \code{disadopt} contract of \code{\link{rdiffnet}}: a function of
#' \code{(expo, cumadopt, time)} that returns a list of length \eqn{Q}
#' (one entry per behaviour) with the integer node indices that
#' disadopt at the current step.
#'
#' Use any of these as the \code{disadopt} argument of \code{rdiffnet},
#' or write your own factory that returns a function with the same
#' signature.
#'
#' @details
#' The four kernels below cover the common cases:
#'
#' \describe{
#'  \item{\code{disadoptmech_random}}{Each currently-adopted node
#'    disadopts independently with probability \code{prob}. Models
#'    constant-rate recovery (the SIR \eqn{\gamma}).}
#'  \item{\code{disadoptmech_bithreshold}}{Currently-adopted nodes
#'    disadopt when their exposure crosses an upper threshold,
#'    \code{threshold_dis}. Pair with \code{adoptmech_threshold} to
#'    instantiate the bi-threshold model of Alipour \emph{et al.}
#'    (2024) — adopt at the lower threshold, disadopt at the upper.}
#'  \item{\code{disadoptmech_logit}}{Bernoulli rule with logit link.
#'    Disadopt with probability
#'    \code{plogis(beta0 + beta_expo * expo)}. To make recovery
#'    \emph{less likely} as exposure grows, set \code{beta_expo < 0}.}
#'  \item{\code{disadoptmech_probit}}{Bernoulli rule with probit link.
#'    Disadopt with probability
#'    \code{pnorm(beta0 + beta_expo * expo)}.}
#' }
#'
#' @param prob Numeric scalar in \eqn{[0, 1]}. Per-step disadoption
#'   probability for \code{disadoptmech_random}.
#' @param threshold_dis Numeric scalar or vector of length \eqn{n}.
#'   Upper-threshold cut-off for \code{disadoptmech_bithreshold}.
#'   A scalar is recycled across all nodes.
#' @param beta0 Numeric scalar. Intercept of the logit/probit
#'   disadoption probability.
#' @param beta_expo Numeric scalar. Slope on exposure for the
#'   logit/probit disadoption probability.
#'
#' @return A function with signature
#'   \code{function(expo, cumadopt, time)} suitable as the
#'   \code{disadopt} argument of \code{\link{rdiffnet}}.
#'
#' @references
#' Alipour, F., Dokshin, F., Maleki, Z., Song, Y., & Ramazi, P. (2024).
#' Enough but not too many: A bi-threshold model for behavioral
#' diffusion. \emph{PNAS Nexus} 3(10).
#' \doi{10.1093/pnasnexus/pgae428}
#'
#' @examples
#' set.seed(2026)
#'
#' # Constant-rate recovery: each adopter recovers with prob 0.10 / step
#' dn <- rdiffnet(n = 50, t = 12, seed.graph = "small-world",
#'                seed.p.adopt = 0.10, stop.no.diff = FALSE,
#'                disadopt = disadoptmech_random(prob = 0.10))
#'
#' # Bi-threshold model (Alipour 2024): adopt when exposure >= 0.30,
#' # disadopt when exposure >= 0.70.
#' dn <- rdiffnet(n = 50, t = 12, seed.graph = "small-world",
#'                seed.p.adopt   = 0.10, stop.no.diff = FALSE,
#'                threshold.dist = 0.30,
#'                disadopt       = disadoptmech_bithreshold(threshold_dis = 0.70))
#'
#' # Logit recovery, exposure-dependent
#' dn <- rdiffnet(n = 50, t = 12, seed.graph = "small-world",
#'                seed.p.adopt = 0.10, stop.no.diff = FALSE,
#'                disadopt = disadoptmech_logit(beta0 = -1, beta_expo = -2))
#'
#' @name disadoption_mechanisms
NULL

#' @rdname disadoption_mechanisms
#' @export
disadoptmech_random <- function(prob) {
  if (missing(prob))
    stop("-disadoptmech_random- requires -prob-.")
  if (!is.numeric(prob) || length(prob) != 1L || is.na(prob) ||
      prob < 0 || prob > 1)
    stop("-prob- must be a single number in [0, 1].")
  force(prob)

  function(expo, cumadopt, time) {
    Q <- dim(cumadopt)[3L]
    out <- vector("list", Q)
    for (q in seq_len(Q)) {
      currently <- which(cumadopt[, time, q] == 1L)
      out[[q]] <- if (length(currently))
        currently[stats::runif(length(currently)) < prob]
      else
        integer()
    }
    out
  }
}

#' @rdname disadoption_mechanisms
#' @export
disadoptmech_bithreshold <- function(threshold_dis) {
  if (missing(threshold_dis))
    stop("-disadoptmech_bithreshold- requires -threshold_dis-.")
  if (!is.numeric(threshold_dis) || anyNA(threshold_dis))
    stop("-threshold_dis- must be a numeric scalar or vector with no NA.")
  force(threshold_dis)

  function(expo, cumadopt, time) {
    n <- dim(cumadopt)[1L]
    Q <- dim(cumadopt)[3L]

    th <- if (length(threshold_dis) == 1L)
      rep(threshold_dis, n)
    else
      threshold_dis
    if (length(th) != n)
      stop("-threshold_dis- length (", length(th),
           ") does not match number of nodes (", n, ").")

    # rdiffnet hands disadopt() an -expo- of shape n x 1 x Q (current
    # time slice only); -cumadopt- carries the full n x T x Q history.
    out <- vector("list", Q)
    for (q in seq_len(Q)) {
      currently <- which(cumadopt[, time, q] == 1L)
      e <- expo[, 1L, q]
      out[[q]] <- currently[e[currently] >= th[currently]]
    }
    out
  }
}

#' @rdname disadoption_mechanisms
#' @export
disadoptmech_logit <- function(beta0, beta_expo) {
  if (missing(beta0) || missing(beta_expo))
    stop("-disadoptmech_logit- requires both -beta0- and -beta_expo-.")
  force(beta0); force(beta_expo)

  function(expo, cumadopt, time) {
    # -expo- shape is n x 1 x Q (current slice only); -cumadopt- is full history.
    Q <- dim(cumadopt)[3L]
    out <- vector("list", Q)
    for (q in seq_len(Q)) {
      currently <- which(cumadopt[, time, q] == 1L)
      if (length(currently)) {
        e <- expo[currently, 1L, q]
        p <- stats::plogis(beta0 + beta_expo * e)
        out[[q]] <- currently[stats::runif(length(p)) < p]
      } else {
        out[[q]] <- integer()
      }
    }
    out
  }
}

#' @rdname disadoption_mechanisms
#' @export
disadoptmech_probit <- function(beta0, beta_expo) {
  if (missing(beta0) || missing(beta_expo))
    stop("-disadoptmech_probit- requires both -beta0- and -beta_expo-.")
  force(beta0); force(beta_expo)

  function(expo, cumadopt, time) {
    # -expo- shape is n x 1 x Q (current slice only); -cumadopt- is full history.
    Q <- dim(cumadopt)[3L]
    out <- vector("list", Q)
    for (q in seq_len(Q)) {
      currently <- which(cumadopt[, time, q] == 1L)
      if (length(currently)) {
        e <- expo[currently, 1L, q]
        p <- stats::pnorm(beta0 + beta_expo * e)
        out[[q]] <- currently[stats::runif(length(p)) < p]
      } else {
        out[[q]] <- integer()
      }
    }
    out
  }
}
