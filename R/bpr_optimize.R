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
#' out_opt <- bpr_optim(x = data, method = "BFGS")
#'
#' @export
bpr_optim <- function(x, ...){
  UseMethod("bpr_optim")
}


# Default function for the generic function 'bpr_optim'
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
#' @return A list containing the following elements:
#' \itemize{
#'  \item{ \code{W_opt}: An Nx(M+1) matrix with the optimized parameter values.
#'    Each row of the matrix corresponds to each element of the list x. The
#'    columns are of the same length as the parameter vector w (i.e. number
#'    of basis functions).
#'  }
#'  \item{ \code{Mus}: An N x M matrix with the RBF centers if basis object is
#'   \code{\link{rbf.object}}, otherwise NULL.}
#'  \item{ \code{basis}: The basis object.
#'  }
#'  \item{ \code{w}: The initial values of the parameters w.
#'  }
#'  \item{ \code{x_extrema}: The min and max values of each promoter region.
#'  }
#' }
#'
#' @seealso \code{\link{bpr_optim}}, \code{\link{bpr_optim.matrix}}
#'
#' @examples
#' data <- bpr_data
#' out_opt <- bpr_optim(x = data, method = "CG")
#'
#' @export
bpr_optim.list <- function(x, w = NULL, basis = NULL, fit_feature = NULL,
                                          method = "CG", itnmax = 100, ...){
  N <- length(x)
  assertthat::assert_that(N > 0)
  out <- do_checks(w = w, basis = basis)
  w   <- out$w
  basis <- out$basis

  # Data frame for storing all the coefficients for each element of list x
  if (is.null(fit_feature)){
    num_features = length(w)
  }else{
    num_features = length(w) + 1
  }
  W_opt <- matrix(NA_real_, nrow = N, ncol = num_features)
  colnames(W_opt) <- paste("w", seq(1, num_features), sep = "")
  # Matrix for storing the centers of RBFs if object class is 'rbf'

  Mus <- NULL
  if (is(basis, "rbf")){
    Mus <- matrix(NA_real_, nrow = N, ncol = basis$M)
    colnames(Mus) <- paste("mu", seq(1, basis$M), sep = "")
  }
  x_extrema <- matrix(NA_real_, nrow = N, ncol = 2)
  # Perform optimization for each element of x, i.e. for each region i.
  for (i in 1:N){
    out_opt <- bpr_optim.matrix(x           = x[[i]],
                                w           = w,
                                basis       = basis,
                                fit_feature = fit_feature,
                                method      = method,
                                itnmax      = itnmax)
    W_opt[i, ] <- out_opt$w_opt
    if (is(basis, "rbf")){
      Mus[i, ] <- out_opt$basis$mus
    }
    x_extrema[i, ] <- out_opt$x_extrema
  }
  return(list(W_opt = W_opt,
              Mus = Mus,
              basis = basis,
              w = w,
              x_extrema = x_extrema))
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
#' @param fit_feature Additional feature on how well the profile fits the
#'  methylation data.
#' @param method The optimization method to be used. See \code{\link[stats]{optim}}
#'  for possible methods. Default is 'CG'.
#' @param itnmax Optional argument giving the maximum number of iterations for
#'  the corresponding method. See \code{\link[stats]{optim}} for details.
#' @param ... Additional parameters
#'
#' @return A list containing the following elements:
#' \itemize{
#'  \item{ \code{w_opt}: Optimized values for the coefficient vector w.
#'    The length of the result is the same as the length of the vector w.
#'  }
#'  \item{ \code{basis}: The basis object.
#'  }
#' }
#'
#' @seealso \code{\link{bpr_optim}}, \code{\link{bpr_optim.list}}
#'
#' @examples
#' basis <- polynomial.object(M=2)
#' w <- c(0.1, 0.1, 0.1)
#' data <- bpr_data[[1]]
#' out_opt <- bpr_optim(x = data, w = w, basis = basis, method = "CG", fit_feature = "NLL")
#'
#' @importFrom stats optim
#'
#' @export
bpr_optim.matrix <- function(x, w = NULL, basis = NULL, fit_feature = NULL,
                                          method = "CG", itnmax = 100, ...){
  obs <- as.vector(x[ ,1])
  data <- x[ ,2:3]
  # Create design matrix H
  des_mat <- design_matrix(x = basis, obs = obs)
  H <- des_mat$H
  basis <- des_mat$basis

  # Call optim function to perform minimization of the NLL of BPR function
  w_opt <- optim(par     = w,
                 fn      = bpr_likelihood,
                 gr      = bpr_derivative,
                 method  = method,
                 control = list(maxit = itnmax),
                 H       = H,
                 data    = data,
                 is_NLL  = TRUE)$par

  # If we need to add the goodness of fit to the data as feature
  if (!is.null(fit_feature)){
    if (identical(fit_feature, "NLL")){
      fit <- bpr_likelihood(w = w_opt,
                            H = H,
                            data = data,
                            is_NLL = TRUE)
    }else if (identical(fit_feature, "RMSE")){
      # Predictions of the target variables
      f_pred <- as.vector(pnorm(H %*% w_opt))
      f_true <- data[,2] / data[, 1]
      fit <- sqrt(mean((f_pred - f_true) ^ 2))
    }
    w_opt <- c(w_opt, fit)
  }

  return(list(w_opt = w_opt,
              basis = basis,
              x_extrema = c(min(obs), max(obs))))
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
    stop("Coefficients vector should be M+1, M: number of basis functions!")
  }
  return(list(w = w, basis = basis))
}
