#' Generic function for performing Gibbs sampling on BPR model
#'
#' \code{bpr_optim} is a generic function which calles the appropriate methods
#' depending on the class of the object \code{x}. Object \code{x} can be either
#' a \code{\link{list}} or a \code{\link{matrix}}.
#'
#' @param x The input object
#' @param ... Additional parameters
#'
#' @seealso \code{\link{bpr_gibbs.list}}, \code{\link{bpr_gibbs.matrix}}
#'
#' @examples
#' data <- bpr_data
#' out_opt <- bpr_gibbs(x = data, is_parallel = FALSE)
#'
#' @export
bpr_gibbs <- function(x, ...){
  UseMethod("bpr_gibbs")
}


# Default function for the generic function 'bpr_gibbs'
bpr_gibbs.default <- function(x, ...){
  stop("Object x should be either matrix or list!")
}


#' Gibbs sampling for the BPR model using list x
#'
#' \code{bpr_gibbs.list} computes the posterior of the BPR model using auxiliary
#' variable approach. Since we cannot compute the posterior analytically, a
#' Gibbs sampling scheme is used. This method calls
#' \code{\link{bpr_gibbs.matrix}} to process each element of the list.
#'
#' @param x A list of elements of length N, where each element is an L x 3
#'  matrix of observations, where 1st column contains the locations. The 2nd
#'  and 3rd columns contain the total trials and number of successes at the
#'  corresponding locations, repsectively.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @inheritParams bpr_gibbs.matrix
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
#' ex_data <- bpr_data
#' basis <- rbf.object(M=3)
#' out_opt <- bpr_gibbs(x = ex_data, is_parallel = FALSE, basis = basis)
#'
#' @export
bpr_gibbs.list <- function(x, w = NULL, basis = NULL, fit_feature = NULL,
                           cpg_dens_feat = FALSE, w_0_mean = NULL,
                           w_0_cov = NULL, gibbs_nsim = 100,
                           is_parallel = TRUE, no_cores = NULL, ...){
  # Check that x is a list object
  assertthat::assert_that(is.list(x))

  # Extract number of observations
  N <- length(x)
  assertthat::assert_that(N > 0)

  # Perform checks for initial parameter values
  out <- .do_checks(w = w, basis = basis)
  w   <- out$w
  basis <- out$basis

  if (is.null(w_0_mean)){
    w_0_mean <- rep(0, length(w))
  }else{
    if (length(w) != length(w_0_mean)){
      w_0_mean <- rep(0, length(w))
    }
  }

  if (is.null(w_0_cov)){
    w_0_cov <- diag(1, length(w))
  }else{
    if (length(w) != NROW(w_0_cov)){
      w_0_cov <- diag(1, length(w))
    }
  }

  # Initialize so the CMD check on R passes without NOTES
  i <- 0

  # If parallel mode is ON
  if (is_parallel){
    # If number of cores is not given
    if (is.null(no_cores)){
      no_cores <- parallel::detectCores() - 2
    }else{
      if (no_cores >= parallel::detectCores()){
        no_cores <- parallel::detectCores() - 1
      }
    }
    if (is.na(no_cores)){
      no_cores <- 2
    }
    # Create cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)

    # Parallel optimization for each element of x, i.e. for each region i.
    res <- foreach::"%dopar%"(obj = foreach::foreach(i = 1:N),
                  ex  = {
                    out_opt <- bpr_gibbs.matrix(x           = x[[i]],
                                                w           = w,
                                                basis       = basis,
                                                fit_feature = fit_feature,
                                                cpg_dens_feat = cpg_dens_feat,
                                                w_0_mean    = w_0_mean,
                                                w_0_cov     = w_0_cov,
                                                gibbs_nsim  = gibbs_nsim)
                              })
    # Stop parallel execution
    parallel::stopCluster(cl)
  }else{
    # Sequential optimization for each element of x, i.e. for each region i.
    res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
                 ex  = {
                   out_opt <- bpr_gibbs.matrix(x           = x[[i]],
                                               w           = w,
                                               basis       = basis,
                                               fit_feature = fit_feature,
                                               cpg_dens_feat = cpg_dens_feat,
                                               w_0_mean    = w_0_mean,
                                               w_0_cov     = w_0_cov,
                                               gibbs_nsim  = gibbs_nsim)
                 })
  }

  # Matrix for storing optimized coefficients
  W_opt <- sapply(res, function(x) x$w_opt)
  if (is.matrix(W_opt)){
    W_opt <- t(W_opt)
  }else{
    W_opt <- as.matrix(W_opt)
  }
  colnames(W_opt) <- paste("w", seq(1, NCOL(W_opt)), sep = "")

  # Matrix for storing the centers of RBFs if object class is 'rbf'
  Mus <- NULL
  if (methods::is(basis, "rbf")){
    if (is.null(basis$mus)){
      Mus <- sapply(lapply(res, function(x) x$basis), function(y) y$mus)
      if (is.matrix(Mus)){
        Mus <- t(Mus)
      }else{
        Mus <- as.matrix(Mus)
      }
      colnames(Mus) <- paste("mu", seq(1, NCOL(Mus)), sep = "")
    }
  }

  # Matrix for storing extrema promoter values
  x_extrema <- t(sapply(res, function(x) x$x_extrema))

  return(list(W_opt = W_opt,
              Mus = Mus,
              basis = basis,
              w = w,
              x_extrema = x_extrema))
}


#' Gibbs sampling for the BPR model using list x
#'
#' \code{bpr_gibbs.matrix} computes the posterior of the BPR model using auxiliary
#' variable approach. Since we cannot compute the posterior analytically, a
#' Gibbs sampling scheme is used.
#'
#' @param x An L x 3 matrix of observations, where 1st column contains the
#'  locations. The 2nd and 3rd columns contain the total trials and number of
#'  successes at the corresponding locations, repsectively.
#' @param w A vector of parameters (i.e. coefficients of the basis functions)
#' @param basis A 'basis' object. See \code{\link{polynomial.object}}
#' @param fit_feature Additional feature on how well the profile fits the
#'  methylation data.
#' @param cpg_dens_feat Additional feature for the CpG density across the
#'  promoter region.
#' @param w_0_mean The prior mean hyperparameter for w
#' @param w_0_cov The prior covariance hyperparameter for w
#' @param gibbs_nsim Optional argument giving the number of simulations of the
#'  Gibbs sampler.
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
#' w_0_mean <- rep(0, length(w))
#' w_0_cov <- diag(10, length(w))
#' data <- bpr_data[[1]]
#' out_opt <- bpr_gibbs(x = data, w = w, w_0_mean = w_0_mean,
#'                      w_0_cov = w_0_cov, basis = basis)
#'
#' basis <- polynomial.object(M=0)
#' w <- c(0.1)
#' w_0_mean <- rep(0, length(w))
#' w_0_cov <- diag(10, length(w))
#' data <- bpr_data[[1]]
#' out_opt <- bpr_gibbs(x = data, w = w, w_0_mean = w_0_mean,
#'                      w_0_cov = w_0_cov, basis = basis)
#'
#' @importFrom stats optim
#' @importFrom truncnorm rtruncnorm
#' @importFrom mvtnorm rmvnorm
#'
#' @export
bpr_gibbs.matrix <- function(x, w = NULL, basis = NULL, fit_feature = NULL,
                             cpg_dens_feat = FALSE, w_0_mean = NULL,
                             w_0_cov = NULL, gibbs_nsim = 100, ...){

  # Vector for storing CpG locations relative to TSS
  obs <- as.vector(x[ ,1])

  # Create design matrix H
  des_mat <- design_matrix(x = basis, obs = obs)
  H       <- des_mat$H
  basis   <- des_mat$basis

  # Total number of reads for each CpG
  N_i     <- as.vector(x[,2])
  # Corresponding number of methylated reads for each CpG
  y_star  <- as.vector(x[,3])

  z <- rep(0, sum(N_i))				   # Latent Normal Variable
  y <- vector(mode = "integer")
  for (i in 1:length(y_star)){
    y <- c(y, rep(1, y_star[i]), rep(0, N_i[i] - y_star[i]))
  }

  N <- sum(N_i)
  n1  <- sum(y)   # Number of successes
  n0  <- N - n1   # Number of failures

  # Create concatenated matrix of observations
  rows <- rep(1:nrow(H), N_i)
  H <- as.matrix(H[rows, ])

  # Store Results
  w_samples <- matrix(0, nrow = gibbs_nsim, ncol = length(w))

  # Start Gibbs sampler
  prec_W_0 <- solve(w_0_cov)
  W_cov <- solve(prec_W_0 + crossprod(H, H))

  for (tt in 2:gibbs_nsim) {
    # Draw Latent Variable, z, from its full conditional, given y
    mus <- H %*% w			# Update Mean of Z

    # Using truncated normal function
    z[y==0] <- rtruncnorm(n = n0, mean = mus[y==0], sd = 1, a = -Inf, b = 0)
    z[y==1] <- rtruncnorm(n = n1, mean = mus[y==1], sd = 1, a = 0, b = Inf)

    # Posterior mean of w|Z
    w_mean <- W_cov %*% (prec_W_0 %*% w_0_mean + crossprod(H, z))
    w_samples[tt, ] <- w <- c(rmvnorm(1, w_mean, as.matrix(W_cov)))
  }

  # Keep summary statistic and get the mean of the Gibbs samples
  if (NCOL(w_samples) == 1){
    w_opt <- mean(w_samples[(ceiling(gibbs_nsim/4)):gibbs_nsim, ])
  }else{
    w_opt <- as.vector(apply(w_samples[(ceiling(gibbs_nsim/4)):gibbs_nsim, ],
                             2, mean))
  }

  # If we need to add the goodness of fit to the data as feature
  if (!is.null(fit_feature)){
    if (identical(fit_feature, "NLL")){
      fit <- bpr_likelihood(w = w_opt,
                            H = H,
                            data = x[ ,2:3],
                            is_NLL = TRUE)
    }else if (identical(fit_feature, "RMSE")){
      # Predictions of the target variables
      f_pred <- as.vector(pnorm(H %*% w_opt))
      f_true <- x[ ,3] / x[ ,2]
      fit <- sqrt(mean( (f_pred - f_true) ^ 2) )
    }
    w_opt <- c(w_opt, fit)
  }

  # Add as feature the CpG density in the promoter region
  if (cpg_dens_feat){
    w_opt <- c(w_opt, length(obs))
  }

  return(list(w_opt = w_opt,
              basis = basis,
              x_extrema = c(min(obs), max(obs))))
}


# Internal function to make all the appropriate type checks.
.do_checks <- function(w = NULL, basis = NULL){
  if (is.null(basis)){
    basis <- rbf.object(M = 3)
  }
  if (is.null(w)){
    w <- rep(0.5, basis$M + 1)
  }
  if (length(w) != (basis$M + 1) ){
    stop("Coefficients vector should be M+1, M: number of basis functions!")
  }
  return(list(w = w, basis = basis))
}
