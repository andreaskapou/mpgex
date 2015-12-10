#' Generic function for creating a desing matrix H
#'
#' This is a generic function which calles the appropriate methods depending
#' on the class of the object \code{x}. Currently only the 'polynomial' class
#' is implemented.
#'
#' @param x A \code{\link{polynomial.object}}
#' @param ... Additional parameters
#'
#' @seealso \code{\link{polynomial_basis}}, \code{\link{rbf_basis}},
#'  \code{\link{design_matrix.polynomial}}, \code{\link{design_matrix.rbf}}
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
#' @return A list containing the design matrix H. The dimensions of the matrix
#' H are N x (M+1), where N is the length of the observations, and M is the
#' degree of the polynomial.
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


#' Creates an RBF design matrix H
#'
#' \code{design_matrix.rbf} creates a design matrix H using radial
#' basis functions of degree M.
#'
#' @param x A 'basis' object.
#' @param obs A vector of observations.
#' @param ... Additional parameters
#'
#' @return A list containing the design matrix \code{H} and the centers of the
#' rbf function \code{mus}. The dimensions of the matrix H are N x (M+1), where
#' N is the length of the observations, and M is the number of radial basis
#' functions. \code{mus} is a vector of length M containg the centers of rbfs.
#'
#' @seealso \code{\link{design_matrix}}, \code{\link{rbf_basis}}
#'
#' @examples
#' obj <- rbf.object(M=2)
#' obs <- c(0,.2,.5)
#' des_mat <- design_matrix(obj, obs)
#'
#' @export
design_matrix.rbf <- function(x, obs, ...){
  assertthat::assert_that(is(x, "rbf"))
  assertthat::assert_that(is.vector(obs))

  ### TODO Evaluate equally spaced mus as an option???

  mus <- NULL
  N   <- length(obs)  # Length of the dataset
  M   <- x$M          # Number of coefficients
  # If number of basis functions is higher than the total observations
  if (M > N - 1){
    stop("Number of basis functions is higher than number of observations!")
  }
  if (M == 0){
    H <- matrix(1, nrow = N, ncol = M + 1)
  }else{
    if (x$eq_spaced_mus){
      if (x$whole_region){
        mus <- seq(-1, 1, length.out = M)
      }else{
        mus <- seq(min(obs), max(obs), length.out = M)
      }
    }else{
      repeat {
        km <- stats::kmeans(obs, M, iter.max = 30, nstart = 10)  # Use K-means
        if (min(km$size) > 0)  # Only accept non-empty clusters
          break
      }
      mus <- km$centers  # RBF centers
    }

    # Convert the 'obs' vector to an N x 1 dimensional matrix
    obs <- as.matrix(obs)
    gamma <- x$gamma  # Inverse width of basis function
    H <- matrix(1, nrow = N, ncol = M)
    for (j in 1:M){
      # TODO: Implement the multivariate case!
      H[ ,j] <- apply(obs, 1, rbf_basis, mus = mus[j], gamma = gamma)
    }
    H <- cbind (1 , H)  # Add the 'bias' term
  }
  return(list(H = H, mus = mus))
}
