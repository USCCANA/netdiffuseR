context("Set of functions used in: plot_diffnet2")

# ------------------------------------------------------------------------------
test_that("round_to_seq", {
  set.seed(12313)

  # Respect the lenght
  x   <- runif(1000, 0,100)
  x_r <- round_to_seq(x, nlevels = 30)
  expect_equal(30L, length(unique(x_r)))
})

# ------------------------------------------------------------------------------
test_that("methods of plot_diffnet2", {
  set.seed(11222)
  dn <- rdiffnet(100, 5)

  # Should get the same outputs
  set.seed(1);ans_dn <- plot_diffnet2(dn)
  set.seed(1);ans_gp <- plot_diffnet2(dn$graph[[nslices(dn)]], dn$toa)
  expect_equal(ans_dn, ans_gp)

})

# ------------------------------------------------------------------------------
test_that("methods of diffMap", {
  set.seed(11445)
  dn <- rdiffnet(50, 5, seed.nodes = "central")

  # Should get the same outputs
  set.seed(1); ans_dn <- diffmap(dn)
  set.seed(1); ans_gp <- diffmap(dn$graph[[nslices(dn)]], dn$toa)
  expect_equal(ans_dn, ans_gp)

  set.seed(123);ans0 <- plot_diffnet2(dn, add.map = "last", layout = ans_gp$coords)
  set.seed(123);ans1 <- diffmap(dn, layout=igraph::norm_coords(ans_gp$coords))

  # Should be the same as adding the map after
  expect_equal(ans0$diffmap, ans1)
  expect_output(print(ans1), "of class.+diffmap")
  expect_output(print(ans1), "of class.+diffmap")
  expect_silent(image(ans1))
  expect_silent(plot(ans1))
})


