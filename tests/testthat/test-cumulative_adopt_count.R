context("Cumulative adopt count, and hazard rate")

toa <- 2001L:2005L

adopt <- toa_mat(toa)
count <- cumulative_adopt_count(adopt$cumadopt)

count_hand <-rbind(
  1:5,
  (1:5)/5
)

count_hand <- as.numeric(rbind(
  count_hand,
  c(0,(count_hand[1,2:5] - count_hand[1,1:4])/count_hand[1,1:4])
))

test_that("Cumulative adopters",
          expect_equal(as.numeric(count), count_hand,
                       tolerance=getOption("diffnet.tol"), check.attributes=FALSE, scale=1)
          )

# Now, hazard rate
hr <- as.numeric(hazard_rate(adopt$cumadopt))
hr_hand <- as.numeric(c(0,(count[1,2:5] - count[1,1:4])/(length(toa) - count[1,1:4])))


test_that("Hazard rate", {
  expect_equal(hr, hr_hand, tolerance=getOption("diffnet.tol"), scale=1)
}
          )
