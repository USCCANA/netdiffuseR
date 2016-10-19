context("Read/write foreign formats")

test_that("Pajek files", {
  # Reading sampson data
  path <- system.file("extdata", "SAMPSON.NET", package = "netdiffuseR")
  SAMPSON <- read_pajek(path)

  g <- edgelist_to_adjmat(SAMPSON$edges[[1]][,1:2], SAMPSON$edges[[1]][,3],
                          recode.ids = FALSE, multiple = FALSE)
  for (i in 2:length(SAMPSON$edges))
    g <- g + edgelist_to_adjmat(SAMPSON$edges[[i]][,1:2],
                                SAMPSON$edges[[i]][,3],
                                recode.ids = FALSE, multiple = FALSE)

  ans0 <- dgr(g, valued = TRUE)
  path <- system.file("extdata", "SAMPSON_W_DEGREE.vec", package="netdiffuseR")
  ans1 <- as.numeric(readLines(path)[-1])

  # Comparing against what pajek says
  expect_equivalent(as.numeric(ans0),ans1)


})

# ------------------------------------------------------------------------------
test_that("DL files", {
  # Reading sampson data
  path <- system.file("extdata", "SAMPSON.DAT", package = "netdiffuseR")
  SAMPSON <- read_ml(path)

  ans0 <- dgr(SAMPSON$adjmat, valued = TRUE)
  ans0 <- rowSums(ans0)

  path <- system.file("extdata", "SAMPSON_W_DEGREE.vec", package="netdiffuseR")
  ans1 <- as.numeric(readLines(path)[-1])

  # Comparing against what pajek says
  expect_equivalent(as.numeric(ans0),ans1)
})

# ------------------------------------------------------------------------------
# test_that("UCINET files", {
#   # Reading sampson data
#   path <- system.file("extdata", "sampson.##d", package = "netdiffuseR")
#   SAMPSON <- read_ucinet(path)ans0 <- dgr(SAMPSON$adjmat, valued = TRUE)
#   ans0 <- rowSums(ans0)
#
#   path <- system.file("extdata", "SAMPSON_W_DEGREE.vec", package="netdiffuseR")
#   ans1 <- as.numeric(readLines(path)[-1])
#
#   # Comparing against what pajek says
#   expect_equivalent(as.numeric(ans0),ans1)
# })
