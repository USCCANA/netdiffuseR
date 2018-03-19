library(netdiffuseR)

data("medInnovationsDiffNet")

diffnet_formula <- function(model) {

  # Get the terms
  ans  <- terms.formula(model)

  # Fetching the diffnet
  obj <- eval(attr(ans, "variables")[[2]])

  if (!inherits(obj, "diffnet"))
    stop("The Left Hand Side of the formula is not a `diffnet` object.",
         call. = FALSE)

  vars <- eval(attr(ans, "variables"))

  # is the first a diffnet object?
  if (!inherits(vars[[1]], "diffnet"))
    stop("Not an object of class diffnet.")

  list(
    1,
    ans
  )
}

ans <- diffnet_formula(medInnovationsDiffNet ~ exposure(valued = TRUE))
ans[[1]][[2]](ans[[1]][[1]])
