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
#' out <- polynomial_basis(X=3, M=2)
#' # Or using a vector X
#' out <- polynomial_basis(X=c(2,3,4), M=2)
#'
#' @export
polynomial_basis <- function(X, M = 1){
  return(X ^ M)
}


#' Apply radial basis function
#'
#' Applied the RBF function to the input X with center mus.
#'
#' @param X Input data.
#' @param mus Centers from where we should compute the distance of the data X.
#' @param gamma Inverse width of radial basis function.
#'
#' @return Input X, after being transformed from the RBF.
#'
#' @examples
#' out <- rbf_basis(X = c(1,2), mus = c(1,1))
#' #
#' out <- rbf_basis(X = c(1,2), mus = c(1,1), gamma = 0.1)
#'
#' @export
rbf_basis <- function(X, mus, gamma = 1){
  return(exp((-1) * gamma * norm(as.matrix(X - mus), "F") ^ 2))
}

##-------------------------------------------------------------

#' Creates a polynomial object
#'
#' \code{polynomial.object} creates objects of type "polynomial". These objects
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
  # Check that M is integer
  assertthat::assert_that(M %% 1 == 0)
  assertthat::assert_that(M > -1)
  obj <- list(M = M)
  class(obj) <- "polynomial"
  return(obj)
}


#' Creates an RBF object
#'
#' \code{rbf.object} creates objects of type "rbf". These objects can be used as
#' input to complex functions in order to perform computation depending on the
#' basis function. E.g see \code{\link{design_matrix}} generic function.
#'
#' @param M The degree of the RBF object that will be created.
#' @param gamma Inverse width of radial basis function.
#'
#' @return An object of type 'RBF'.
#'
#' @seealso \code{\link{design_matrix}}
#'
#' @examples
#' obj <- rbf.object(M=2)
#' (obj)
#'
#' @export
rbf.object <- function(M = 2, gamma = 1){
  # Check that M is integer
  assertthat::assert_that(M %% 1 == 0)
  assertthat::assert_that(M > -1)
  obj <- list(M = M, gamma = gamma)
  class(obj) <- "rbf"
  return(obj)
}
