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

  n <- 50
  num_of_behaviors <- 1

  # Must work

  x <- 0.35 # numeric scalar
  thr <- rdiffnet_make_threshold(x, n = n, num_of_behaviors = num_of_behaviors)
  expect_equal(dim(thr), c(n, 1))

  x <- runif(n) # vector of length n
  thr <- rdiffnet_make_threshold(x, n = n, num_of_behaviors = num_of_behaviors)
  expect_equal(dim(thr), c(n, 1))

  x <- function() runif(1) # function
  thr <- rdiffnet_make_threshold(x, n = n, num_of_behaviors = num_of_behaviors)
  expect_equal(dim(thr), c(n, 1))


  # Must show ERROR

  x <- runif(100)# Length greater than n
  expect_error(
    rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors),
    "Incorrect threshold input in function -rdiffnet_make_threshold-."
  )

  x <- runif(25)# Length less than n
  expect_error(
    rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors),
    "Incorrect threshold input in function -rdiffnet_make_threshold-."
  )

  x <- "invalid_input"# Non-numeric input
  expect_error(
    rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors),
    "Incorrect threshold input in function -rdiffnet_make_threshold-."
  )

})

# Multiple --------------------------------------------------------------------

test_that("Multi diff models rdiff args work", {

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

test_that("Checking threshold for multiple diffusion", {

  n <- 50
  num_of_behaviors <- 2

  # Must work

  x <- matrix(runif(100), nrow = n, ncol = num_of_behaviors) # matrix input
  thr <- rdiffnet_make_threshold(x, n = n, num_of_behaviors = num_of_behaviors)
  expect_equal(dim(thr), c(n, num_of_behaviors))

  x <- list(function() runif(1), function() rexp(1)) # list of functions
  thr <- rdiffnet_make_threshold(x, n = n, num_of_behaviors = num_of_behaviors)
  expect_equal(dim(thr), c(n, num_of_behaviors))

  x <- list(0.14,0.05) # list of scalars
  thr <- rdiffnet_make_threshold(x, n = n, num_of_behaviors = num_of_behaviors)
  expect_equal(dim(thr), c(n,num_of_behaviors))

  x <- list(runif(n), runif(n)) # list of vectors
  thr <- rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors)
  expect_equal(dim(thr),c(n,num_of_behaviors))


  # Must show ERROR

  x <- list(runif(2*n),runif(n)) # incorrect vector length (too long for one behavior)
  expect_error(
    rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors),
    "Incorrect threshold input in function -rdiffnet_make_threshold-."
  )

  x <- list(runif(n/2),runif(n)) # incorrect vector length (too short for one behavior)
  expect_error(
    rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors),
    "Incorrect threshold input in function -rdiffnet_make_threshold-."
  )

  x <- c(runif(n),runif(n)) # the input should be a list
  expect_error(
    rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors)
  )

  x <- list(0.14) # Only one behavior provided in the list
  expect_error(
    rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors),
    "The length of the list must match the number of behaviors"
  )

  x <- runif(n) # Only one behavior provided in the vector
  expect_error(
    rdiffnet_make_threshold(x,n=n,num_of_behaviors=num_of_behaviors)
  )
})
