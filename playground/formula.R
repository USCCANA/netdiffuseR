library(netdiffuseR)

data("medInnovationsDiffNet")

diffreg <- function(model) {

  # Checks ---------------------------------------------------------------------

  # It must be a formula
  if (!inherits(model, "formula"))
    stop("Not a formula. `model` should be a formula in the form of <diffnet object> ~ exposure + ...", call. = FALSE)

  # Get the terms and checking whether exposure is present
  terms <- terms.formula(model)
  labs  <- attr(terms, "term.labels")
  test  <- which(grepl("^exposure(\\(|$)", labs))

  if (!length(test))
    stop("No `exposure` term in the formula. The `diffreg` should include `exposure` on the RHS of the formula.", call. = FALSE)

  exposure_term <- labs[test]


  # Fetching the diffnet
  LHS <- attr(terms, "variables")[[2]]

  if (!inherits(eval(LHS), "diffnet"))
    stop("The LHS of the formula is not a `diffnet` object. `diffreg` is only for `diffnet` objects.",
         call. = FALSE)

  # Updating the formula -------------------------------------------------------
  # Now the model will be Adopt ~ Exposure + ...
  model_call <- update.formula(model, paste("~", exposure_term, "+."))
  model <- update.formula(
    update.formula(model, paste("Adopt ~ + . -", exposure_term)),
    ~ exposure + .
    )


  # Computing exposure ---------------------------------------------------------
  exposure_term <- attr(terms, "variables")[[test + 2L]]

  # Replacing the network object
  if (!is.call(exposure_term))
    exposure_term <- as.call(list(exposure_term, graph = LHS))
  else
    exposure_term$graph <- LHS

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

  eval(
    call(
      "glm", noquote(model), data=quote(dat),
      family = quote(binomial(link = "logit")),
      subset = quote(ifelse(is.na(toa), TRUE, toa >= per))
      )
    )

}

ans <- diffreg(
  medInnovationsDiffNet ~ factor(city) + proage + exposure(lags = 1L) + per)
summary(ans)

ans <- diffreg(
  kfamilyDiffNet ~ pregs + exposure(lags = 1L) + per)
summary(ans)


ans <- diffreg(
  brfarmersDiffNet ~ age + I(age^2) + exposure(lags = 1L) + factor(per) + income)
summary(ans)

