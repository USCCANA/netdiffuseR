context("Infection and susceptibility")

################################################################################
# Susceptibility tests
################################################################################

test_that("Numeric tests of Susceptibility", {
  # Preparing data -------------------------------------------------------------
  graph <- array(0,dim=c(5,5,4), dimnames = list(1:5,1:5,1:4))
  graph[5,1:4,] <- 1
  times <- c(1, 4, 1, 2, 3)

  graphl <- lapply(1:4, function(x) methods::as(graph[,,x], "dgCMatrix"))
  names(graphl) <- 1:4

  gl <- list(`Array`=graph, `List`=graphl)

  for(i in names(gl)) {
    # K=1 ------------------------------------------------------------------------
    suscep <- netdiffuseR::susceptibility(gl[[i]], times, normalize = FALSE)[5,1]
    expect_equal(suscep, 1/3, scale=1, tol=getOption("diffnet.tol"), info=i)

    # K=2 ------------------------------------------------------------------------
    suscep <- netdiffuseR::susceptibility(gl[[i]], times, K = 2, normalize=FALSE)[5,1]
    expect_equal(suscep, 2/4, scale=1, tol=getOption("diffnet.tol"), info=i)

    # K=1 exponential r=.1--------------------------------------------------------
    suscep <- netdiffuseR::susceptibility(gl[[i]], times, normalize = FALSE,
                                          r=.1, expdiscount = TRUE)[5,1]
    expect_equal(suscep, 1/3, scale=1, tol=getOption("diffnet.tol"), info=i)

    # K=2 exponential r=.1--------------------------------------------------------
    bmark <- (1/1 + 2/(1.1))/(3/1 + 2/(1.1))
    suscep <- netdiffuseR::susceptibility(gl[[i]], times, K = 2, normalize=FALSE,
                                          r=.1, expdiscount = TRUE)[5,1]
    expect_equal(suscep, bmark, scale=1, tol=getOption("diffnet.tol"), info=i)
  }
}
)



################################################################################
# Infection tests
################################################################################
test_that("Numeric tests of Susceptibility", {
  # Preparing data -------------------------------------------------------------
  graph <- array(0,dim=c(5,5,4), dimnames = list(1:5,1:5,1:4))
  graph[2:5,1,] <- 1
  times <- c(1, 4, 2, 3, 3) +1
  times[2] <- 1

  graphl <- lapply(1:4, function(x) methods::as(graph[,,x], "dgCMatrix"))
  names(graphl) <- 1:4

  gl <- list(`Array`=graph, `List`=graphl)

  for(i in names(gl)) {
    # K=1 ------------------------------------------------------------------------
    infect <- netdiffuseR::infection(gl[[i]], times, normalize = FALSE)
    expect_equal(infect[1,1], 1/3, scale=1, tol=getOption("diffnet.tol"), info=i)

    # K=2 ------------------------------------------------------------------------
    infect <- netdiffuseR::infection(gl[[i]], times, normalize = FALSE, K=2)
    expect_equal(infect[1,1], 1/2, scale=1, tol=getOption("diffnet.tol"), info=i)

    # K=1 exponential r=.1--------------------------------------------------------
    infect <- netdiffuseR::infection(gl[[i]], times, normalize = FALSE,
                                          r=.1, expdiscount = TRUE)
    expect_equal(infect[1,1], 1/3, scale=1, tol=getOption("diffnet.tol"), info=i)

    # K=2 exponential r=.1--------------------------------------------------------
    bmark <- (1/1 + 2/(1.1))/(3/1 + 2/(1.1))
    infect <- netdiffuseR::infection(gl[[i]], times, K = 2, normalize=FALSE,
                                          r=.1, expdiscount = TRUE)
    expect_equal(infect[1,1], bmark, scale=1, tol=getOption("diffnet.tol"), info=i)
  }
}
)
