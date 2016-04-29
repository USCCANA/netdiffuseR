context("Foreign function")

test_that("Igraph back and forth", {
  # Getting the data
  g <- medInnovationsDiffNet
  ig <- diffnet_to_igraph(g)
  dn <- igraph_to_diffnet(ig[[1]], "toa")

  g$graph  <- as.array(g)
  dn$graph <- as.array(dn)
  expect_equal(g, dn)
})
