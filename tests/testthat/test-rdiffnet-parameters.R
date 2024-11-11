# Single --------------------------------------------------------------------

# Must work
test_that(
  "Checking single diffusion rdiffnet args", {

  # Must work

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

test_that("Checking threshold for single diffusion", {

  # Must work

  thr <- rdiffnet_make_threshold(1.5, n = 50, q = 1)
  expect_equal(dim(thr), c(50, 1))

  x <- runif(50)
  thr <- rdiffnet_make_threshold(x, n = 50, q = 1)
  expect_equal(dim(thr), c(50, 1))

  thr <- rdiffnet_make_threshold(function() 0.5, n = 50, q = 1)
  expect_equal(dim(thr), c(50, 1))

  thr <- rdiffnet_make_threshold(function() rexp(1), n = 50, q = 1)
  expect_equal(dim(thr), c(50, 1))

  # Must show ERROR

  x <- runif(100) # Length n*q
  expect_error(
    rdiffnet_make_threshold(x, n = 50, q = 1)
  )

})

# Multiple --------------------------------------------------------------------

test_that("Multi diff  models rdiff args work", {

  # Must work

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


# NOT working now !!!

# test_that("Checking threshold for multiple diffusion", {
#
#   # Must work
#
# x <- matrix(runif(100), nrow = 50, ncol = 2)
# thr <- rdiffnet_make_threshold(x, n = 50, q = 2)
# expect_equal(dim(thr), c(50, 2))
#
# x <- runif(100) # Length n*q
# expect_error(
#   rdiffnet_make_threshold(x, n = 50, q = 2)
# )
#
#   seed.p.adopt <- list(function() runif(1), function() rexp(1))
#   thr <- rdiffnet_make_threshold(seed.p.adopt, n = 50, q = 2)
#   expect_equal(dim(thr), c(50,1))
#
#
#   seed.p.adopt <- list(0.14,0.05)
#   thr <- rdiffnet_make_threshold(seed.p.adopt[[1]], n = 50, q = 2)
#   expect_equal(dim(thr), c(50,2))
#
#
#   seed.p.adopt <- list(runif(50), runif(50))
#
#   # Test first element of list
#   thr <- rdiffnet_make_threshold(seed.p.adopt[[1]], n = 50, q =1 )
#
#   expect_equal(dim(thr), c(50,1))
#
#
#   # Must show ERROR
#
#   x <- runif(100) # Length n*q
#   expect_error(
#     rdiffnet_make_threshold(x, n=100,q=3),
#     "incorrect input
# }
