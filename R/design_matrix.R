#' Generic function for creating a desing matrix H
#'
#' This is a generic function which calles the appropriate methods depending
#' on the class of the object \code{x}. Currently only the 'polynomial' class
#' is implemented.
#'
#' @param x A \code{\link{polynomial.object}}
#' @param ... Additional parameters
#'
#' @seealso \code{\link{polynomial_basis}}, \code{\link{design_matrix.polynomial}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' obs <- c(0,.2,.5)
#' H <- design_matrix(obj, obs)
#'
#' @export
design_matrix <- function(x, ...){
  UseMethod("design_matrix")
}


# Default function for the generic function 'design_matrix'
design_matrix.default <- function(x, ...){
  stop(paste("Object of type '", class(x), "' is not implemented.", sep = ""))
}


#' Creates a polynomial design matrix H
#'
#' \code{design_matrix.polynomial} creates a design matrix H using polynomial
#' basis functions of degree M.
#'
#' @param x A 'basis' object.
#' @param obs A vector of observations.
#' @param ... Additional parameters
#'
#' @return The design matrix H. The dimensions of the matrix H are
#' N x (M+1), where N is the length of the observations, and M is the degree of
#' the polynomial.
#'
#' @seealso \code{\link{design_matrix}}, \code{\link{polynomial_basis}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' obs <- c(0,.2,.5)
#' H <- design_matrix(obj, obs)
#'
#' @export
design_matrix.polynomial <- function(x, obs, ...){
  assertthat::assert_that(is(x, "polynomial"))
  assertthat::assert_that(is.vector(obs))

  N   <- length(obs)  # Length of the dataset
  M   <- x$M + 1      # Number of coefficients
  H <- matrix(1, nrow = N, ncol = M)
  for (j in 1:M){
    H[ ,j] <- polynomial_basis(obs, j - 1)  # Compute X^(j-1)
  }
  return(list(H = H))
}
