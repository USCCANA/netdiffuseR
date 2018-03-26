context("Diffusion Regression")

test_that("Errors and warnings", {

  data("medInnovationsDiffNet")
  expect_error(diffreg(~exposure), "not a `diffnet`")
  expect_error(diffreg(medInnovationsDiffNet~proage), "No `exposure`")

  expect_warning(
    diffreg(medInnovationsDiffNet~proage+exposure(lags = 0) + factor(per)),
    "not lagged"
  )
  expect_warning(
    diffreg(medInnovationsDiffNet~proage+exposure),
    "naturally increases over time"
  )


})
