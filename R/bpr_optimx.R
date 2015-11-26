#' Generic function for optimizing the BPR log likelihood function
#'
#' This is a generic function which calles the appropriate methods depending
#' on the class of the object \code{x}. Currently only the 'polynomial' class
#' is implemented.
#'
#' @param x The input object
#' @param ... Additional parameters
#'
#' @seealso \code{\link{bpr_optimx.list}}, \code{\link{bpr_optimx.matrix}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' obs <- c(0,.2,.5)
#' H <- design_matrix(obj, obs)
#'
#' @export
bpr_optimx <- function(x, ...){
  UseMethod("bpr_optimx")
}


bpr_optimx.default <- function(x, ...){
  stop("Object x should be either matrix or list!")
}

#' Optimization method for the BPR log likelihood function using a list x
#'
#' \code{bpr_optimx.list} computes the maximum log likelihood of the BPR
#' function. Since it cannot be evaluated analytically, an optimization
#' procedure is used. The \code{\link[optimx]{optimx}} packages is used for
#' performing optimization. This method calls \code{\link{bpr_optimx.matrix}}
#' to process each element of the list.
#'
#' @param x A list of elements of length N, where each element is an L x 3
#'  matrix of observations, where 1st column contains the locations. The 2nd
#'  and 3rd columns contain the total trials and number of successes at the
#'  corresponding locations, repsectively.
#' @inheritParams bpr_optimx.matrix
#'
#' @return An N x M matrix, with the optimized parameter values. Each row of
#'  the matrix corresponds to each element of the list x. The columns are of
#'  the same length as the parameter vector w (i.e. number of basis functions).
#'
#' @seealso \code{\link{bpr_optimx}}, \code{\link{bpr_optimx.matrix}}
#'
#' @export
bpr_optimx.list <- function(x, w = NULL, basis = NULL, method="CG", itnmax = 200, ...){
  N <- length(x)
  assertthat::assert_that(N > 0)

  out <- make_checks(w = w, basis = basis)
  w   <- out$w
  basis <- out$basis

  # Matrix for storing all the coefficients for each element of list x
  W <- matrix(0, nrow = N, ncol = length(w))

  # Perform optimization for each element of list x,
  # i.e. for each promoter region i.
  for (i in 1:N){
    W[i, ] <- bpr_optimx.matrix(x = x[[i]],
                                w = w,
                                basis = basis,
                                method = method,
                                itnmax = itnmax)
  }
  return(W)
}

#' Optimization method for the BPR log likelihood function using a matrix x
#'
#' \code{bpr_optimx.matrix} computes the maximum log likelihood of the BPR
#' function. Since it cannot be evaluated analytically, an optimization
#' procedure is used. The \code{\link[optimx]{optimx}} packages is used for
#' performing optimization.
#'
#' @param x An L x 3 matrix of observations, where 1st column contains the
#'  locations. The 2nd and 3rd columns contain the total trials and number of
#'  successes at the corresponding locations, repsectively.
#' @param w A vector of parameters (i.e. coefficients of the basis functions)
#' @param basis A 'basis' object. See \code{\link{polynomial.object}}
#' @param method The optimization method to be used. See \code{\link[optimx]{optimx}}
#'  for possible methods. Default is 'CG'.
#' @param itnmax Optional argument giving the maximum number of iterations for
#'  the corresponding method. See \code{\link[optimx]{optimx}} for details.
#' @param ... Additional parameters
#'
#' @return Optimized values for the parameter/coefficient vector w. The length
#'  of the result is the same as the length of the vector w.
#'
#' @seealso \code{\link{bpr_optimx}}, \code{\link{bpr_optimx.list}}
#'
#'
#' @export
bpr_optimx.matrix <- function(x, w = NULL, basis = NULL, method="CG", itnmax = 200, ...){
  # Make all the appropriate type checks
  out <- make_checks(w = w, basis = basis)
  w <- out$w
  basis <- out$basis

  obs <- as.vector(x[ ,1])
  data <- x[ ,2:3]
  # Create design matrix H
  des_mat <- design_matrix(x = basis, obs = obs)
  H <- des_mat$H

  # Call optimx function to perform optimization
  opt_w <- optimx::optimx(par     = w,
                          fn      = bpr_likelihood,
                          gr      = bpr_derivative,
                          method  = method,
                          itnmax  = itnmax,
                          control = list(maximize = TRUE),
                          H       = H,
                          data    = data)

  return(opt_w[1:length(w)])
}


# Internal function to make all the appropriate type checks.
make_checks <- function(w, basis){
  if (is.null(basis)){
    basis <- polynomial.object()
  }
  if (is.null(w)){
    w <- rep(0.1, basis$M + 1)
  }
  if (length(w) != (basis$M + 1) ){
    stop("Length of coefficients vector should be M + 1, where M is the number of basis functions!")
  }
  return(list(w = w, basis = basis))
}
