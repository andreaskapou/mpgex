context("Checking utils functions work")

test_that("minmax scaling works",{
  data <- c(-20, 0, 15, 20)
  xmin <- -20
  xmax <- 20
  expect_identical(minmax_scaling(data, xmin, xmax), c(0, 0.5, 0.875, 1))
  data <- c(-20, 20)
  expect_identical(minmax_scaling(data, fmin=-2, fmax=2), c(-2, 2))
  expect_identical(minmax_scaling(data, -xmin, xmax), c(-20, 20))
  expect_identical(minmax_scaling(4, 0, 8), 0.5)
})
