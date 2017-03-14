#' Matching Estimators with Network Data
#'
#' \strong{WARNING}: This function is still in development and has not been tested throughly.
#' Following Aral et al. (2009), \code{netmatch} computes matching
#' estimators for network data. The function \code{netmatch_prepare}, which
#' prepares the data to be used with \code{\link[MatchIt:matchit]{matchit}} from
#' the \pkg{\link[MatchIt:MatchIt]{MatchIt}} package, is called by \code{netmatch}.
#'
#' @param dat \code{data.frame} with dynamic data. Must be of
#' \code{nrow(dat)==nslices(graph)*nnodes(graph)}.
#' @param graph List with sparse matrices.
#' @param timevar Character scalar. Name of time variable
#' @param depvar Character scalar. Name of the dependent variable
#' @param covariates Character vector. Name(s) of the control variable(s).
#' @param adopt_thr Either a numeric scalar or vector of length \code{nslices(graph)}.
#' Sets the threshold of \code{depvar} at which it is considered that an observation
#' has adopted a behavior.
#' @param treat_thr Either a numeric scalar or vector of length \code{nslices(graph)}.
#' Sets the threshold of \code{exposure} at which it is considered that an
#' observation is treated.
#' @param expo_pcent Logical scalar. When \code{TRUE}, exposure is computed
#' non-normalized (so it is a count rather than a percentage).
#' @param expo_lag Integer scalar. Number of lags to consider when computing
#' exposure. \code{expo_lag=1} defines exposure in T considering behavior and
#' network at T-1.
#' @param ... Further arguments to be passed to \code{\link[MatchIt:matchit]{matchit}}.
#' @return In the case of \code{netmatch_prepare}
#' \item{dat}{A \code{data.frame} with the original data (covariates), plus the
#' following new variables: \code{treat}, \code{adopt}, \code{exposure}.
#' }
#' \item{match_model}{A formula to be passed to \code{netmatch}}
#'
#' \code{netmatch} returns the following:
#' \item{fATT}{A numeric vector of length \eqn{N_1}{N1} (number of treated used
#' in the matching process). Treatment effects on the treated at the individual
#' level}
#' \item{match_obj}{The output from \code{matchit}.}
#' @name netmatch
#' @author
#' George G. Vega Yon
#' @references
#' Aral, S., Muchnik, L., & Sundararajan, A. (2009). Distinguishing
#' influence-based contagion from homophily-driven diffusion in dynamic networks.
#' Proceedings of the National Academy of Sciences of the United States of America,
#' 106(51), 21544–21549. \doi{10.1073/pnas.0908800106}
#'
#' Imbens, G. W., & Wooldridge, J. M. (2009). Recent Developments in the
#' Econometrics of Program Evaluation. Journal of Economic Literature, 47(1),
#' 5–86. \doi{10.1257/jel.47.1.5}
#'
#' King, G., & Nielsen, R. (2015). Why Propensity Scores Should Not Be Used for.
#'
#' Sekhon, J. S. (2008). The Neyman-Rubin Model of Causal Inference and Estimation
#' Via Matching Methods. The Oxford Handbook of Political Methodology.
#' \doi{10.1093/oxfordhb/9780199286546.003.0011}
#' @details
#'
#' In Aral et al. (2009), the matching estimator is used as a response to the
#' fact that the observed network is homophilous. Essentially, using exposure
#' as a treatment indicator, which is known to be endogenous, we can apply the
#' same principle of matching estimators in which, after controlling for characteristics
#' (covariates), individuals from the treated group (exposed to some behavior)
#' can be compared to individuals from the control group (not exposed to that
#' behavior), as the only difference between the two is the exposure.
#'
#' As pointed out in King & Nielsen (2015), it is suggested that, contrary to
#' what Aral et al. (2009), the matching is not performed over propensity score
#' since it is know that the later can increase imbalances in the data and thus
#' obtaining exactly the opposed outcome that matching based estimators pursue.
#'
#' A couple of good references for matching estimators are Imbens and Wooldridge
#' (2009), and Sekhon (2008).
NULL

#' @export
#' @rdname netmatch
netmatch_prepare <- function (
  dat, graph, timevar, depvar, covariates,
  treat_thr = rep(1L, length(graph)),
  adopt_thr = rep(1L, length(graph)),
  expo_pcent = FALSE,
  expo_lag = 0L
) {

  # Subset
  Y    <- as.matrix(dat[, depvar, drop=FALSE])
  X    <- as.matrix(dat[, covariates, drop=FALSE])
  Time <- as.matrix(dat[, timevar, drop=FALSE])

  # Checking thresholds
  if (length(treat_thr)==1)
    treat_thr <- rep(treat_thr, length(graph))
  if (length(adopt_thr)==1)
    adopt_thr <- rep(adopt_thr, length(graph))

  # Obtaining constants
  pers  <- sort(unique(Time))
  npers <- length(pers)
  n     <- length(Y)/npers

  if (expo_lag >= npers)
    stop("Not enought time points for such lag.")

  # Generating exposures and treatment
  expo  <- matrix(NA, ncol=npers, nrow=n)
  adopt <- matrix(NA, ncol=npers, nrow=n)
  treat <- matrix(NA, ncol=npers, nrow=n)
  ans   <- vector("list", npers)
  for (i in seq_along(pers)) {

    # Subset data
    per <- pers[i]
    y   <- Y[Time == per,, drop=FALSE]
    x   <- X[Time == per,, drop=FALSE]

    # Computing exposure and seting treatment
    adopt[,i] <- as.integer(y >= adopt_thr[i])

    # Skiping lags
    if (i <= expo_lag) next

    # Computing exposure (and normalizing if required)
    expo[,i]  <- as.vector(graph[[i-expo_lag]] %*% adopt[,i-expo_lag,drop=FALSE])
    if (expo_pcent)
      expo[,i] <- expo[,i]/(Matrix::rowSums(graph[[i-expo_lag]]) + 1e-20)

    # Assigning regime
    treat[,i] <- as.integer(expo[,i] >= treat_thr[i])

    # Creating dataset
    ans[[i]] <- data.frame(Y=y, X=x, treat=treat[,i],
                           adopt = adopt[,i-expo_lag],
                           expo  = expo[,i], check.names = FALSE)
    colnames(ans[[i]]) <-
      c("Y", covariates, "treat", "adopt", "expo")

  }

  # Model
  match_model <- as.formula(paste("treat ~", paste0(covariates, collapse=" + ")))

  # Returning outcome
  return(list(dat=ans, match_model=match_model))
}


#' @rdname netmatch
#' @export
netmatch <- function(
  dat, graph, timevar, depvar, covariates,
  treat_thr  = rep(1L, length(graph)),
  adopt_thr  = rep(1L, length(graph)),
  expo_pcent = FALSE,
  expo_lag   = 0L, ...) {

  # No so good for now.
  if (nslices(graph) != 2)
    stop("-netmatch- is currently supported only for nslices(graph) == 2.")

  # Preparing the data
  ans   <- netmatch_prepare(dat=dat, graph = graph, timevar = timevar,
                            depvar = depvar, covariates = covariates,
                            treat_thr = treat_thr, adopt_thr = adopt_thr,
                            expo_pcent = expo_pcent, expo_lag = expo_lag)

  # print(str(ans[[1]]))
  # Processing the data so we can use it with -matchit-
  dat <- ans$dat[[2]]
  rownames(dat) <- 1:nrow(dat)

  # Matching (with replacement)
  match_obj   <- MatchIt::matchit(ans$match_model, data=dat, ...)

  mmethod <- match_obj$call$method
  if (length(mmethod) && ("cem" %in% mmethod)) {

    # Retrieving sub classes
    dat1 <- MatchIt::match.data(match_obj, "treat")
    dat0 <- MatchIt::match.data(match_obj, "control")

    # Computing effect. At the end CEM gives weights (very simple)
    # coef(lm(Y~treat, data=match.data(match_obj), weights=weights))["treat"]
    fATT <- with(dat1, sum(Y*weights)/sum(weights)) -
      with(dat0, sum(Y*weights)/sum(weights))

  } else {
    match_index <- match_obj$match.matrix

    # Computing feasible Average Treatment Effect on the Treated
    fATT <- dat[,"Y"]

    matched_control <- lapply(1:nrow(match_index), function(x) {
      x <- match_index[x,,drop=TRUE]
      as.integer(x[!is.na(x)])
    })

    # Computing differences
    fATT <-Map(function(i1, i0) mean(fATT[i1] - fATT[i0]),
               i0=as.integer(rownames(match_index)),
               i1=matched_control
    )

    fATT <- mean(unlist(fATT, recursive = TRUE), na.rm=TRUE)
  }


  # Wrapping up
  return(list(
    fATT      = fATT,
    match_obj = match_obj,
    exposures = dat$expo
  )
  )

}
