#' Diffusion regression model
#'
#' A wrapper of \code{glm}, this function estimates a lagged regression model of
#' adoption as a function of exposure and other controls as especified by the
#' user.
#'
#' @param model An object of class formula where the right-hand-side is an object of
#' class \code{\link[=diffnet]{diffnet}}
#' @param type Character scalar. Either \code{"probit"} or \code{"logit"}.
#' @details
#' The model must be in the following form:
#'
#' \preformatted{
#' <diffnet object> ~ exposure + covariate1 + covariate2 + ...
#' }
#'
#' Where \code{exposure} can be especified either as a simple term, or as a
#' call to the exposure function, e.g. to compute exposure with a lag of
#' length 2, the formula could be:
#'
#' \preformatted{
#' <diffnet object> ~ exposure(lags = 2) + covariate1 + covariate2 + ...
#' }
#'
#' When no argument is passed to \code{exposure}, the function sets a lag
#' of length 1 by default (see the \emph{Lagged regression} section).
#'
#' This is a wrapper of \code{\link[stats:glm]{glm}}. The function does the
#' following steps:
#' \enumerate{
#'  \item Compute exposure by calling \code{exposure} on the LHS (dependent variable).
#'  \item Modify the formula so that the model is on adoption as a function of
#'  exposure and whatever covariates the user specifies.
#'  \item Selects either \code{"probit"} or \code{"logit"} and prepares the call
#'  to \code{glm}. This includes passing the following line:
#'  \preformatted{
#'  subset = ifelse(is.na(toa), TRUE, toa >= per)
#'  }
#'  This results in including observations that either did not adopted or up to
#'  the time of adoption.
#'  \item Estimates the model.
#' }
#'
#' The data passed to \code{glm} is obtained by using \code{\link{as.data.frame.diffnet}}.
#'
#' @section Lagged regression:
#'
#' The model estimated is a lagged regression model that has two main assumptions:
#' \enumerate{
#' \item The network is exogenous to the behavior (no selection effect)
#' \item The influence effect (diffusion) happens in a lagged fasion, hence,
#' exposure is computed lagged.
#' }
#'
#' If either of these two assumptions is not met, then the model becomes endogenous,
#' ans so inference becomes invalid.
#'
#' In the case of the first assumption, the user can overcome the non-exogeneity
#' problem by providing an alternative network. This can be done by especifying
#' \code{alt.graph} in the \code{exposure} function so that the network becomes
#' exogenous to the adoption.
#'
#' @return
#' An object of class \code{\link[stats:glm]{glm}}.
#' @examples
#' data("medInnovationsDiffNet")
#'
#' # Default model
#' ans <- diffreg(
#'   medInnovationsDiffNet ~ exposure + factor(city) + proage + per)
#' summary(ans)
#' @export
diffreg <- function(model, type=c("logit", "probit")) {

  # Checks ---------------------------------------------------------------------

  # It must be a formula
  if (!inherits(model, "formula"))
    stop("Not a formula. `model` should be a formula in the form of <diffnet object> ~ exposure + ...", call. = FALSE)

  # Get the terms and checking whether exposure is present
  terms <- stats::terms.formula(model)
  labs  <- attr(terms, "term.labels")
  test  <- which(grepl("^exposure(\\(|$)", labs))

  if (!length(test))
    stop("No `exposure` term in the formula. The `diffreg` should include `exposure` on the RHS of the formula.", call. = FALSE)

  exposure_term <- labs[test]


  # Fetching the diffnet
  LHS <- attr(terms, "variables")[[2]]

  if (!inherits(eval(LHS), "diffnet"))
    stop("The LHS of the formula (dependent variable) is not a `diffnet` object. `diffreg` is only for `diffnet` objects.",
         call. = FALSE)

  # Updating the formula -------------------------------------------------------
  # Now the model will be Adopt ~ Exposure + ...
  model_call <- stats::update.formula(model, paste("~", exposure_term, "+."))
  model <- stats::update.formula(
    stats::update.formula(model, paste("Adopt ~ + . -", exposure_term)),
    ~ exposure + .
    )


  # Computing exposure ---------------------------------------------------------
  exposure_term <- attr(terms, "variables")[[test + 2L]]

  # Replacing the network object
  if (!is.call(exposure_term)) {
    exposure_term <- as.call(list(exposure_term, graph = LHS, lags = 1L))
  } else {
    if (is.null(exposure_term$lags) || exposure_term$lags < 1L)
      warning(
        "When the model is not lagged, i.e. is not `exposure(..., lags = 1L)`, ",
        "the model becomes endogenous, so inference is invalid."
      )
    exposure_term$graph <- LHS
  }


  exposure_term <- eval(exposure_term)

  dat <- eval(LHS)
  dat[["exposure"]] <- exposure_term
  dat[["Adopt"]]    <- dat$cumadopt
  dat <- as.data.frame(dat)

  # Some warnings
  if (!any(grepl("^(factor\\()?per", labs)))
    warning(
      "The variable `per` is not included in the model. ",
      "This can bias the `exposure` as adoption naturally increases over time.",
      call. = FALSE
      )

  # Making the call ------------------------------------------------------------

  ans <- call(
    "glm", noquote(model), data=quote(dat),
    subset = quote(ifelse(is.na(toa), TRUE, toa >= per))
    )

  if (all(type == c("logit", "probit")))
    type <- "logit"

  if (type == "probit")
    ans$family <- quote(binomial(link = "probit"))
  else if (type == "logit")
    ans$family <- quote(binomial(link = "logit"))
  else
    stop("The model of type `", type, "` is not supported.", call. = FALSE)

  eval(ans)

}
