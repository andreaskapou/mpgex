#' Apply polynomial basis function.
#'
#' Applies the polynomial basis function of degree M to the input X.
#'
#' @param X The input data, either a scalar, vector or matrix.
#' @param M Integer, denoting the degree of the polynomial basis that will be
#'  applied to the input X. M should not be negative.
#'
#' @return Input X, after being transformed from the polynomial basis function.
#'
#' @examples
#' data <- polynomial_basis(X=3, M=2)
#' # Or using a vector X
#' data <- polynomial_basis(X=c(2,3,4), M=2)
#'
#' @export
polynomial_basis <- function(X, M){
  return(X ^ M)
}


#' Creates a polynomial object
#'
#' \code{polynomial.basis} creates objects of type "polynomial". These objects
#' can be used as input to complex functions in order to perform computation
#' depending on the basis function. E.g see \code{\link{design_matrix}} generic
#' function.
#'
#' @param M The degree of the polynomial object that will be created.
#'
#' @return An object of type 'polynomial'
#'
#' @seealso \code{\link{design_matrix}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' (obj)
#'
#' @export
polynomial.object <- function(M = 1){
  assertthat::assert_that(M %% 1 == 0)
  assertthat::assert_that(M > -1)
  obj <- list(M = M)
  class(obj) <- "polynomial"
  return(obj)
}
