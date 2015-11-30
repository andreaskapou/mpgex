context("Performing calculations on Basis functions")

test_that("polynomial function works fine", {
  expect_identical(polynomial_basis(c(4), 2), 16)
  expect_identical(polynomial_basis(c(2,6), 2), c(4, 36))
  expect_identical(polynomial_basis(matrix(c(1,2,3,4), ncol = 2), 2),
                   matrix(c(1, 4, 9, 16), ncol = 2))
})

test_that("polynomial basis object works fine", {
  expect_error(polynomial.object(1.5))
  expect_error(polynomial.object(-1))
  expect_is(polynomial.object(2), "polynomial")
})

test_that("design matrix works fine", {
  polyn <- polynomial.object(M = 2)
  obs <- c(1,2,3)
  expect_is(design_matrix(polyn, obs), "list")
  expect_error(design_matrix(polyn, cbind(1, obs)))
  expect_identical(design_matrix(polyn, obs),
                   list(H=matrix(c(1,1,1,1,2,3,1,4,9), ncol=3)))

  wr_basis <- list(M=2)
  class(wr_basis) <- "basis_x"
  expect_error(design_matrix(wr_basis))
})
