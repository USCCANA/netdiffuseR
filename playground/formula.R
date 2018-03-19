library(netdiffuseR)

data("medInnovationsDiffNet")

diffnet_formula <- function(model) {

  # Checks ---------------------------------------------------------------------

  # It must be a formula
  if (!inherits(model, "formula"))
    stop("Not a formula", call. = FALSE)

  # Get the terms and checking whether exposure is present
  terms <- terms.formula(model)
  labs  <- attr(terms, "term.labels")
  test  <- which(grepl("^exposure(\\(|\\+|\\s)", labs))

  if (!length(test))
    stop("No `exposure` term in the formula", call. = FALSE)

  exposure_term <- labs[test]


  # Fetching the diffnet
  LHS <- attr(terms, "variables")[[2]]

  if (!inherits(eval(LHS), "diffnet"))
    stop("The Left Hand Side of the formula is not a `diffnet` object.",
         call. = FALSE)

  diffnet_term <- labs[1L]

  # Updating the formula -------------------------------------------------------
  # Now the model will be Adopt ~ Exposure + ...
  model_call <- update.formula(model, paste("~", exposure_term, "+."))
  model <- update.formula(model, paste("Adopt ~ exposure + . -", exposure_term))

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

  # Filtering the data ---------------------------------------------------------

  glm(model, data=dat, family = gaussian(link = "logit"))

}

ans <- diffnet_formula(medInnovationsDiffNet ~ factor(city) + proage + exposure(valued = TRUE))
