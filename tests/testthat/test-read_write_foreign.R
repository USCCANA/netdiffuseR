context("Read/write foreign formats")

test_that("Pajek files", {
  # Reading sampson data
  path <- system.file("extdata", "SAMPSON.NET", package = "netdiffuseR")
  SAMPSON <- read_pajek(path)

  g <- edgelist_to_adjmat(SAMPSON$edges[[1]][,1:2],
                          recode.ids = FALSE, multiple = TRUE)
  for (i in 2:length(SAMPSON$edges))
    g <- g + edgelist_to_adjmat(SAMPSON$edges[[i]][,1:2],
                                recode.ids = FALSE, multiple = TRUE)
})
