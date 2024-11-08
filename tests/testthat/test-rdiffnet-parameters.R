# Single --------------------------------------------------------------------

# Must work
test_that(
  "Checking single diffusion rdiffnet args", {
  seed.p.adopt <- c(0.14)
  seed.nodes <- c('random')
  behavior <- c("random behavior")
  rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)

  class(rdiffnet_args$seed.p.adopt) == "list"
  class(rdiffnet_args$seed.nodes) == "list"
  class(rdiffnet_args$behavior) == "list"

  seed.p.adopt <- 0.14
  seed.nodes <- 'random'
  behavior <- "random behavior"
  rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  rdiffnet_args$seed.p.adopt; class(rdiffnet_args$seed.p.adopt)
  rdiffnet_args$seed.nodes; class(rdiffnet_args$seed.nodes)
  rdiffnet_args$behavior; class(rdiffnet_args$behavior)

  class(rdiffnet_args$seed.p.adopt) == "list"
  class(rdiffnet_args$seed.nodes) == "list"
  class(rdiffnet_args$behavior) == "list"

  # Must show ERROR

  seed.p.adopt <- c(0.4,0.82)
  seed.nodes <- c('random')
  behavior <- "random behavior"
  expect_error(
    rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  )

  seed.p.adopt <- c(0.14)
  seed.nodes <- c('random', 'central')
  behavior <- "random behavior"
  expect_error(
    rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  )

  seed.p.adopt <- c(0.14)
  seed.nodes <- "central"
  behavior <- c("random behavior_1", "random behavior_2")
  expect_error(
    rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  )
  behavior <- list("random behavior_1", "random behavior_2")
  expect_error(
    rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  )
})


# Multiple --------------------------------------------------------------------

# Must work
test_that("Multi diff  models rdiff args work", {
  seed.p.adopt <- list(0.14,0.05)
  seed.nodes <- list('random', "central")
  behavior <- list("random behavior_1", "random behavior_2")
  rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  class(rdiffnet_args$seed.p.adopt) == "list"
  class(rdiffnet_args$seed.nodes) == "list"
  class(rdiffnet_args$behavior) == "list"

  behavior <- c("random behavior_1", "random behavior_2")
  rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  class(rdiffnet_args$behavior) == "list"

  behavior <- "random behavior" #Default
  rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  class(rdiffnet_args$behavior) == "list"

  seed.nodes <- c(1,3,5)
  behavior <- list("random behavior_1", "random behavior_2")
  rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  class(rdiffnet_args$seed.nodes) == 'list'

  seed.nodes <- list('marginal',"central")
  rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  class(rdiffnet_args$seed.nodes) == 'list'

  seed.nodes <- list('marginal',"central") ######
  behavior <- c("random behavior_1")
  rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)

  # Must show ERROR

  seed.p.adopt <- c(0.14,0.05)
  seed.nodes <- list('random', "central")
  behavior <- list("random behavior_1", "random behavior_2")
  expect_error(
    rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  )

  seed.p.adopt <- list(0.14,0.05)
  seed.nodes <- c('marginal',"central")
  expect_error(
    rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  )



  behavior <- list("random behavior_1")
  expect_error(
    rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  )

  seed.nodes <- list('random')
  behavior <- list("random behavior_1", "random behavior_2")
  expect_error(
    rdiffnet_args <- rdiffnet_validate_args(seed.p.adopt, seed.nodes, behavior)
  )
})

