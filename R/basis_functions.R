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
#' @param eq_spaced_mus Logical, if TRUE, equally spaced centers are created,
#'  otherwise centers are created using \code{\link[stats]{kmeans}} algorithm
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
rbf.object <- function(M = 2, gamma = 1, eq_spaced_mus = FALSE){
  # Check that M is integer
  assertthat::assert_that(M %% 1 == 0)
  assertthat::assert_that(M > -1)
  obj <- list(M = M, gamma = gamma, eq_spaced_mus = eq_spaced_mus)
  class(obj) <- "rbf"
  return(obj)
}


#------------------------------------------------------------


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
#' @seealso \code{\link{rbf_basis}}
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
#' @seealso \code{\link{polynomial_basis}}
#'
#' @examples
#' out <- rbf_basis(X = c(1,2), mus = c(1,1))
#' #
#' out <- rbf_basis(X = c(1,2), mus = c(1,1), gamma = 0.1)
#'
#' @export
rbf_basis <- function(X, mus, gamma = 1){
  return(exp( (-1) * gamma * norm(as.matrix(X - mus), "F") ^ 2))
}


#-------------------------------------------------------------


#' Evaluate probit polynomial function
#'
#' Method for evaluating the probit transformation of the polynomial function
#' of degree M for observation data x and coefficients w.
#'
#' @inheritParams eval_polyn_func
#'
#' @return The probit transformed polynomial function values.
#'
#' @examples
#' x <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_prob_polyn_func(x, w)
#'
#' @export
eval_prob_polyn_func <- function(x, w){
  assertthat::assert_that(is.vector(x))
  assertthat::assert_that(is.vector(w))

  f <- eval_polyn_func(x, w)
  Phi <- pnorm(f)
  return(Phi)
}


#' Evaluate probit rbf function
#'
#' Method for evaluating the probit transformation of the rbf function
#' with M basis for observation data x and coefficients w.
#'
#' @inheritParams eval_polyn_func
#' @param basis The basis object.
#' @param mus The centers of the rbf function.
#'
#' @return The probit transformed rbf function values.
#'
#' @examples
#' x <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' basis <- rbf.object(M = 2)
#' mus <- c(2,2.5)
#' out <- eval_prob_rbf_func(x, w, basis, mus)
#'
#' @export
eval_prob_rbf_func <- function(x, w, basis, mus){
  assertthat::assert_that(is.vector(w))

  f <- eval_rbf_func(x, w, basis, mus)
  Phi <- pnorm(f)
  return(Phi)
}


#--------------------------------------------------------


#' Evaluate polynomial function
#'
#' Method for evaluating the the polynomial function of degree M for
#' observation data x and coefficients w.
#'
#' @param x Input / observation data.
#' @param w Vector of length M, containing the coefficients of an Mth-order
#'  basis function.
#'
#' @return The polynomial function values.
#'
#' @examples
#' x <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' out <- eval_polyn_func(x, w)
#'
#' @export
eval_polyn_func <- function(x, w){
  f <- rep(0, length(x))
  M <- length(w)
  for (i in 1:M){
    f <- f + w[i] * polynomial_basis(x, i - 1)
  }
  return(f)
}


#' Evaluate rbf function
#'
#' Method for evaluating the rbf function with M basis for observation
#' data x and coefficients w.
#'
#' @inheritParams eval_polyn_func
#' @param basis The basis object.
#' @param mus The centers of the rbf function.
#'
#' @return The rbf function values.
#'
#' @examples
#' x <- c(1,2,3)
#' w <- c(0.1, 0.3, -0.6)
#' basis <- rbf.object(M = 2)
#' mus <- c(2,2.5)
#' out <- eval_rbf_func(x, w, basis, mus)
#'
#' @export
eval_rbf_func <- function(x, w, basis, mus){
  f <- rep(w[1], length(x))
  M <- basis$M
  x <- as.matrix(x)
  for (i in 1:M){
    f <- f + w[i + 1] * apply(x, 1, rbf_basis,
                              mus = mus[i],
                              gamma = basis$gamma)
  }
  return(f)
}
