#' Generic function for optimizing BPR negative log likelihood function
#'
#' \code{bpr_optim} is a generic function which calles the appropriate methods
#' depending on the class of the object \code{x}. Object \code{x} can be either
#' a \code{\link{list}} or a \code{\link{matrix}}.
#'
#' @param x The input object
#' @param ... Additional parameters
#'
#' @seealso \code{\link{bpr_optim.list}}, \code{\link{bpr_optim.matrix}}
#'
#' @examples
#' data <- bpr_data
#' W_opt <- bpr_optim(x = data, method = "BFGS")
#'
#' @export
bpr_optim <- function(x, ...){
  UseMethod("bpr_optim")
}


bpr_optim.default <- function(x, ...){
  stop("Object x should be either matrix or list!")
}

#' Optimization method for the BPR NLL function using list x
#'
#' \code{bpr_optim.list} minimizes the negative log likelihood of the BPR
#' function. Since it cannot be evaluated analytically, an optimization
#' procedure is used. The \code{\link[stats]{optim}} packages is used for
#' performing optimization. This method calls \code{\link{bpr_optim.matrix}}
#' to process each element of the list.
#'
#' @param x A list of elements of length N, where each element is an L x 3
#'  matrix of observations, where 1st column contains the locations. The 2nd
#'  and 3rd columns contain the total trials and number of successes at the
#'  corresponding locations, repsectively.
#' @inheritParams bpr_optim.matrix
#'
#' @return An N x M matrix, with the optimized parameter values. Each row of
#'  the matrix corresponds to each element of the list x. The columns are of
#'  the same length as the parameter vector w (i.e. number of basis functions).
#'
#' @seealso \code{\link{bpr_optim}}, \code{\link{bpr_optim.matrix}}
#'
#' @examples
#' data <- bpr_data
#' W_opt <- bpr_optim(x = data, method = "CG")
#'
#' @export
bpr_optim.list <- function(x, w = NULL, basis = NULL, method = "CG", itnmax = 100, ...){
  N <- length(x)
  assertthat::assert_that(N > 0)

  out <- do_checks(w = w, basis = basis)
  w   <- out$w
  basis <- out$basis

  # Data frame for storing all the coefficients for each element of list x
  W <- matrix(NA_real_, nrow = N, ncol = length(w))
  colnames(W) <- paste("w", seq(1,length(w)), sep = "")

  # Perform optimization for each element of list x,
  # i.e. for each promoter region i.
  for (i in 1:N){
    W[i, ] <- bpr_optim.matrix(x = x[[i]],
                               w = w,
                               basis = basis,
                               method = method,
                               itnmax = itnmax)
  }
  return(W)
}

#' Optimization method for the BPR NLL using matrix x
#'
#' \code{bpr_optim.matrix} minimizes the negative log likelihood of the BPR
#' function. Since it cannot be evaluated analytically, an optimization
#' procedure is used. The \code{\link[stats]{optim}} packages is used for
#' performing optimization.
#'
#' @param x An L x 3 matrix of observations, where 1st column contains the
#'  locations. The 2nd and 3rd columns contain the total trials and number of
#'  successes at the corresponding locations, repsectively.
#' @param w A vector of parameters (i.e. coefficients of the basis functions)
#' @param basis A 'basis' object. See \code{\link{polynomial.object}}
#' @param method The optimization method to be used. See \code{\link[stats]{optim}}
#'  for possible methods. Default is 'CG'.
#' @param itnmax Optional argument giving the maximum number of iterations for
#'  the corresponding method. See \code{\link[stats]{optim}} for details.
#' @param ... Additional parameters
#'
#' @return Optimized values for the parameter/coefficient vector w. The length
#'  of the result is the same as the length of the vector w.
#'
#' @seealso \code{\link{bpr_optim}}, \code{\link{bpr_optim.list}}
#'
#' @examples
#' basis <- polynomial.object(M=2)
#' w <- c(0.1, 0.1, 0.1)
#' data <- bpr_data[[1]]
#' w_opt <- bpr_optim(x = data, w = w, basis = basis, method = "CG")
#'
#' @importFrom stats optim
#'
#' @export
bpr_optim.matrix <- function(x, w = NULL, basis = NULL, method = "CG", itnmax = 100, ...){

  obs <- as.vector(x[ ,1])
  data <- x[ ,2:3]
  # Create design matrix H
  des_mat <- design_matrix(x = basis, obs = obs)
  H <- des_mat$H

  # Call optim function to perform minimization of the NLL of BPR function
  opt_w <- optim(par     = w,
                 fn      = bpr_likelihood,
                 gr      = bpr_derivative,
                 method  = method,
                 control = list(maxit = itnmax),
                 H       = H,
                 data    = data,
                 is_NLL  = TRUE)$par

  return(as.matrix(opt_w))
}


# Internal function to make all the appropriate type checks.
do_checks <- function(w, basis){
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