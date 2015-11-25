#' Optimization method for the BPR log likelihood function
#'
#'

bpr_optimx <- function(X, w = NULL, basis = NULL, method="CG", itnmax = 200){
  UseMethod("bpr_optimx")
}


bpr_optimx.default <- function(X, ...){
  stop("Object X should be either matrix or list!")
}

#' Optimization method for the BPR log likelihood function using a list X
#'
bpr_optimx.list <- function(X, w = NULL, basis = NULL, method="CG", itnmax = 200){
  N <- length(X)
  assertthat::assert_that(N > 0)

  out <- make_checks(w = w, basis = basis)
  w <- out$w
  basis <- out$basis

  # Matrix for storing all the coefficients for each element of list X
  W <- matrix(0, nrow = N, ncol = length(w))

  # Perform optimization for each element of list X,
  # i.e. for each promoter region i.
  for (i in 1:N){
    W[i, ] <- bpr_optimx.matrix(X = X[[i]],
                                w = w,
                                basis = basis,
                                method = method,
                                itnmax = itnmax)
  }
  return(W)
}


bpr_optimx.matrix <- function(X, w = NULL, basis = NULL, method="CG", itnmax = 200){
  # Make all the appropriate type checks
  out <- make_checks(w = w, basis = basis)
  w <- out$w
  basis <- out$basis

  obs <- as.vector(X[ ,1])
  data <- X[ ,2:3]
  # Create design matrix H
  des_mat <- design_matrix(x = basis, obs = obs)
  H <- des_mat$H

  # Call optimx function to perform optimization
  opt_w <- optimx(par     = w,
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
