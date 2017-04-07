context("Bass Diffusion model")

set.seed(1122333)
dat <- rdiffnet(200, 20)

# ------------------------------------------------------------------------------
test_that("bass fits the right thing", {
  ans0 <- coef(fitbass(dat))

  dat <- cumulative_adopt_count(dat)["prop",]
  Times <- 1:length(dat)

  # Objective function
  objf <- function(par, Times, dat) {
    sum((dat - bass_F(Times, par[1], par[2]))^2)
  }

  # Gradient of objf
  grad <- function(par, Times, dat) {
    g   <- bass_dF(par[1], par[2], Times)
    val <- (bass_F(Times, par[1], par[2]) - dat)
    c(sum(g[1,]*val), sum(g[2,]*val))
  }

  # Fitting using optim
  ans1 <- stats::optim(
    c(p = dat[1], q = .5),
    objf, grad, method="L-BFGS-B", lower = 1e-5, Times = Times, dat = dat
    )$par

  # Expectations
  expect_equal(unname(ans1), unname(ans0), tol = .001)
})

test_that("diffnet_bass methods", {
  ans <- fitbass(dat)

  expect_silent(plot(ans))
  expect_output(print(ans), "Nonlinear")
})
