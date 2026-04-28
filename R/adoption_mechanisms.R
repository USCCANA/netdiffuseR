#' Adoption mechanisms for \code{rdiffnet}
#'
#' A family of pluggable kernels that decide which nodes adopt at each
#' simulation step. Pass any of these as the \code{adoption_mechanism}
#' argument of \code{\link{rdiffnet}}, or write your own function that
#' follows the same contract.
#'
#' @param expo Numeric vector of length \eqn{n}. Per-node exposure at
#'   the current time step (\code{expo[, , q]} for behaviour \eqn{q}).
#' @param thresholds Numeric vector of length \eqn{n}. Per-node adoption
#'   threshold (\code{thr[, q]}). Used by the deterministic kernel;
#'   passed but ignored by the stochastic kernels so user-defined
#'   mechanisms can choose whether to use it.
#' @param not_adopted Logical vector of length \eqn{n}. \code{TRUE} for
#'   nodes that have not yet adopted behaviour \eqn{q}
#'   (\code{is.na(toa[, q])}).
#' @param time Integer scalar. Current simulation time step.
#' @param pars Named list of mechanism-specific parameters. Each
#'   kernel documents which fields it expects.
#'
#' @return Integer vector of node indices that adopt at this step.
#'
#' @details
#' The contract is intentionally minimal so that any user can write a
#' mechanism without reading the package internals: receive the current
#' state, decide who adopts, return the indices. The three kernels
#' below cover the common cases.
#'
#' \describe{
#'  \item{\code{adoptmech_threshold}}{Tom Valente's deterministic
#'    threshold rule. Adopt iff \code{expo[i] >= thresholds[i]}.
#'    Ignores \code{pars}.}
#'  \item{\code{adoptmech_logit}}{Bernoulli rule with logit link.
#'    Adopt with probability \code{plogis(beta0 + beta_expo * expo[i])}.
#'    Requires \code{pars$beta0} and \code{pars$beta_expo}.}
#'  \item{\code{adoptmech_probit}}{Bernoulli rule with probit link.
#'    Adopt with probability \code{pnorm(beta0 + beta_expo * expo[i])}.
#'    Requires \code{pars$beta0} and \code{pars$beta_expo}.}
#' }
#'
#' @examples
#' set.seed(2026)
#'
#' # Default deterministic threshold
#' dn <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
#'                seed.p.adopt = 0.1, stop.no.diff = FALSE)
#'
#' # Stochastic logit mechanism
#' dn <- rdiffnet(n = 30, t = 6, seed.graph = "small-world",
#'                seed.p.adopt = 0.1, stop.no.diff = FALSE,
#'                adoption_mechanism = adoptmech_logit,
#'                adoption_pars      = list(beta0 = -2, beta_expo = 5))
#'
#' @name adoption_mechanisms
NULL

#' @rdname adoption_mechanisms
#' @export
adoptmech_threshold <- function(expo, thresholds, not_adopted, time, pars) {
  which((expo >= thresholds) & not_adopted)
}

#' @rdname adoption_mechanisms
#' @export
adoptmech_logit <- function(expo, thresholds, not_adopted, time, pars) {
  if (is.null(pars$beta0) || is.null(pars$beta_expo))
    stop("-adoptmech_logit- requires -pars- with both -beta0- and -beta_expo-.")
  p <- stats::plogis(pars$beta0 + pars$beta_expo * expo)
  which((stats::runif(length(p)) < p) & not_adopted)
}

#' @rdname adoption_mechanisms
#' @export
adoptmech_probit <- function(expo, thresholds, not_adopted, time, pars) {
  if (is.null(pars$beta0) || is.null(pars$beta_expo))
    stop("-adoptmech_probit- requires -pars- with both -beta0- and -beta_expo-.")
  p <- stats::pnorm(pars$beta0 + pars$beta_expo * expo)
  which((stats::runif(length(p)) < p) & not_adopted)
}
