context("Cumulative adopt count, and hazard rate")

# Generating the data ----------------------------------------------------------

# Common data
set.seed(4912)
toa <- sample(c(2001L:2005L, NA), 10, TRUE)
nper <- max(toa, na.rm = TRUE) - min(toa, na.rm = TRUE) + 1
graph <- lapply(1:nper, function(x) rgraph_ba(t=9))

# Creating a diffnet
diffnet <- as_diffnet(graph, toa, undirected = TRUE)

# Calculating numbers
adopt <- toa_mat(toa)
count <- cumulative_adopt_count(adopt$cumadopt)
count_dn <- cumulative_adopt_count(diffnet)

# Cumadopt ---------------------------------------------------------------------
test_that("Cumulative adopters", {

  # Manual calculations
  count_hand <- sapply(diffnet$meta$pers, function(x) sum(!is.na(toa) & toa <= x))
  count_hand <-rbind(
    count_hand,
    count_hand/length(toa)
  )

  count_hand <- as.numeric(rbind(
    count_hand,
    c(0,(count_hand[1,2:nper] - count_hand[1,1:(nper-1)])/count_hand[1,1:(nper-1)])
  ))

  expect_equal(
    as.numeric(count), count_hand, tolerance=getOption("diffnet.tol"),
    check.attributes=FALSE, scale=1, info="Using default")
  expect_equal(
    as.numeric(count_dn), count_hand, tolerance=getOption("diffnet.tol"),
    check.attributes=FALSE, scale=1, info="Using diffnet object")

})

# Hazard rate ------------------------------------------------------------------
test_that("Hazard rate", {

  hr <- as.numeric(hazard_rate(adopt$cumadopt))
  hr_dn <- as.numeric(hazard_rate(diffnet))

  hr_hand <- as.numeric(
    c(0,(count[1,2:nper] - count[1,1:(nper-1)])/(length(toa) - count[1,1:(nper-1)])))

  expect_equal(hr, hr_hand, tolerance=getOption("diffnet.tol"), scale=1)
  expect_equal(hr_dn, hr_hand, tolerance=getOption("diffnet.tol"), scale=1)
})

# Hazard rate -- status-aware (M11) -------------------------------------------
# Bit-identical for absorbing inputs (legacy guarantee) and SIR-correct for
# non-monotone status arrays (multi-cycle / recovery).

test_that("hazard_rate on a status array with recovery counts re-adoptions", {
  # Hand-built status with two distinct adoption episodes for node 1.
  #   t = 1 2 3 4 5 6
  #   n1: 1 1 0 0 1 1   (adopt t=1, recover t=3, re-adopt t=5)
  #   n2: 0 1 1 0 0 0   (adopt t=2, recover t=4)
  #   n3: 0 0 0 0 0 0   (never adopts)
  #   n4: 0 0 1 1 1 1   (adopt t=3, absorbing)
  status <- rbind(
    c(1L, 1L, 0L, 0L, 1L, 1L),
    c(0L, 1L, 1L, 0L, 0L, 0L),
    c(0L, 0L, 0L, 0L, 0L, 0L),
    c(0L, 0L, 1L, 1L, 1L, 1L)
  )

  # Hand-computed expected hazard:
  #   t=2: fresh = n2 (0->1) = 1; susc at t=1 = n2,n3,n4 = 3; haz = 1/3
  #   t=3: fresh = n4 (0->1) = 1; susc at t=2 = n3,n4 = 2; haz = 1/2
  #   t=4: fresh = 0;             susc at t=3 = n2,n3   = 2; haz = 0
  #   t=5: fresh = n1 (0->1 re-adopt) = 1; susc at t=4 = n1,n2,n3 = 3; haz = 1/3
  #   t=6: fresh = 0;             susc at t=5 = n2,n3   = 2; haz = 0
  expected <- c(0, 1/3, 1/2, 0, 1/3, 0)

  hr <- as.numeric(hazard_rate(status, no.plot = TRUE))
  expect_equal(hr, expected, tolerance = getOption("diffnet.tol"), scale = 1)
})

test_that("hazard_rate via diffnet stays bit-identical for absorbing inputs", {
  # The legacy tests above already cover the absorbing case end-to-end against
  # a hand-built formula. This extra guard confirms the two code paths
  # (diffnet vs raw matrix) produce numerically identical results on the same
  # underlying status.
  hr_via_dn  <- as.numeric(hazard_rate(diffnet, no.plot = TRUE))
  hr_via_mat <- as.numeric(hazard_rate(diffnet$status, no.plot = TRUE))
  expect_equal(hr_via_dn, hr_via_mat,
               tolerance = getOption("diffnet.tol"), scale = 1)
})
