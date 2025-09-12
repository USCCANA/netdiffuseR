context("Testing package .onLoad and options")
test_that("There should be diffnet options", {
  netdiffuseR:::.onLoad()
  expect_is(getOption("diffnet.self"), "logical")
  expect_is(getOption("diffnet.undirected"), "logical")
  expect_is(getOption("diffnet.multiple"), "logical")
  expect_is(getOption("diffnet.tol"), "numeric")
})
