#' Generic function for creating a design matrix H
#'
#' This is a generic function which calles the appropriate methods depending
#' on the class of the object \code{x}.
#'
#' @param x A basis function object.
#' @param ... Additional parameters.
#'
#' @seealso \code{\link{polynomial_basis}}, \code{\link{rbf_basis}},
#'  \code{\link{design_matrix.polynomial}}, \code{\link{design_matrix.rbf}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' obs <- c(0,.2,.5)
#' polyn <- design_matrix(obj, obs)
#'
#' #----------------
#'
#' obj <- rbf.object(M=2)
#' obs <- c(0,.2,.5)
#' rbf <- design_matrix(obj, obs)
#'
#' @export
design_matrix <- function(x, ...){
  UseMethod("design_matrix")
}


# Default function for the generic function 'design_matrix'
design_matrix.default <- function(x, ...){
  stop(paste("Object of type '", class(x), "' is not implemented.", sep = ""))
}


#' Create polynomial design matrix H
#'
#' \code{design_matrix.polynomial} creates a design matrix H using polynomial
#' basis functions of degree M.
#'
#' @param x A basis object.
#' @param obs A vector of observations.
#' @param ... Additional parameters
#'
#' @return A list containing the design matrix H and the basis object. The
#' dimensions of the matrix H are N x (M+1), where N is the length of the
#' observations, and M is the degree of the polynomial.
#'
#' @seealso \code{\link{design_matrix}}, \code{\link{polynomial_basis}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' obs <- c(0,.2,.5)
#' des_mat <- design_matrix(obj, obs)
#'
#' @export
design_matrix.polynomial <- function(x, obs, ...){
  assertthat::assert_that(is(x, "polynomial"))
  assertthat::assert_that(is.vector(obs))

  N <- length(obs)  # Length of the dataset
  H <- matrix(1, nrow = N, ncol = x$M + 1)
  if (x$M > 0){
    for (j in 1:x$M){
      H[ ,j + 1] <- polynomial_basis(obs, j)  # Compute X^(j)
    }
  }
  return(list(H = H, basis = x))
}


#' Creates an RBF design matrix H
#'
#' \code{design_matrix.rbf} creates a design matrix H using radial
#' basis functions of degree M.
#'
#' @inheritParams design_matrix.polynomial
#'
#' @return A list containing the design matrix \code{H} and the basis object.
#' The dimensions of the matrix H are Nx(M+1), where N is the length of the
#' observations, and M is the number of radial basis functions. The updated
#' \code{basis} object contains also the updated centers of RBFs.
#'
#' @seealso \code{\link{design_matrix}}, \code{\link{rbf_basis}}
#'
#' @examples
#' obj <- rbf.object(M=3)
#' obs <- c(0,.2,.5, 0.3)
#' des_mat <- design_matrix(obj, obs)
#'
#' @export
design_matrix.rbf <- function(x, obs, ...){
  assertthat::assert_that(is(x, "rbf"))
  assertthat::assert_that(is.vector(obs))

  N   <- length(obs)  # Length of the dataset
  if (x$M > N - 1){
    stop("Number of basis functions should be less than observations!")
  }
  if (x$M == 0){
    H <- matrix(1, nrow = N, ncol = 1)
    x$mus <- 0
    x$gamma <- 1
  }else{
    if (is.null(x$mus)){
      if (x$eq_spaced_mus){
        x$mus <- vector(mode = "numeric", x$M)
        if (x$whole_region){
          # TODO: Should this be deleted or get in another way the min, max?
          for (i in 1:x$M){
            x$mus[i] <- i * (1 - (-1)) / (x$M + 1) + (-1)
          }
        }else{
          for (i in 1:x$M){
            x$mus[i] <- i * (max(obs) - min(obs)) / (x$M + 1) + min(obs)
          }
        }
      }else{
        repeat{
          km <- stats::kmeans(obs, x$M, iter.max = 30, nstart = 10)
          if (min(km$size) > 0){
            break  # Only accept non-empty clusters
          }
        }
        x$mus <- km$centers  # RBF centers
      }
    }
    if (is.null(x$gamma)){
      if (x$whole_region){
        x$gamma <- x$M^2 / (abs(1) + abs(-1)) ^ 2
      }else{
        x$gamma <- x$M^2 / (abs(max(obs)) + abs(min(obs))) ^ 2
      }
    }
    # Convert the 'obs' vector to an N x 1 dimensional matrix
    obs <- as.matrix(obs)
    H <- matrix(1, nrow = N, ncol = x$M + 1)
    for (j in 1:x$M){
      H[ ,j + 1] <- apply(obs, 1, rbf_basis, mus = x$mus[j], gamma = x$gamma)
    }
  }
  return(list(H = H, basis = x))
}
