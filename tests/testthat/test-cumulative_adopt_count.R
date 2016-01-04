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
