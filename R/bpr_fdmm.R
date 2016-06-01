#' Gibbs sampling algorithm for BPR mixture model
#'
#' \code{bpr_fdmm} implements the Gibbs sampling algorithm for performing
#'  clustering on DNA methylation profiles, where the observation model is the
#'  Binomial distributed Probit Regression function,
#'  \code{\link{bpr_likelihood}}.
#'
#' @param x A list of elements of length N, where each element is an L x 3
#' matrix of observations, where 1st column contains the locations. The 2nd
#' and 3rd columns contain the total trials and number of successes at the
#' corresponding locations, repsectively.
#' @param K Integer denoting the number of clusters K.
#' @param pi_k Vector of length K, denoting the mixing proportions.
#' @param w A MxK matrix, where each column contains the basis function
#' coefficients for the corresponding cluster.
#' @param basis A 'basis' object. E.g. see \code{\link{polynomial.object}}
#' @param w_0_mean The prior mean hyperparameter for w
#' @param w_0_cov The prior covariance hyperparameter for w
#' @param dir_a The Dirichlet concentration parameter, prior over pi_k
#' @param gibbs_nsim Argument giving the number of simulations of the
#'  Gibbs sampler.
#' @param gibbs_burn_in Argument giving the burn in period of the
#'  Gibbs sampler.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#' @param is_verbose Logical, print results during EM iterations
#'
#' @importFrom stats rmultinom
#' @importFrom MCMCpack rdirichlet
#' @export
bpr_fdmm <- function(x, K = 2, pi_k = NULL, w = NULL, basis = NULL,
                     w_0_mean = NULL, w_0_cov = NULL, dir_a = NULL,
                     gibbs_nsim = 5000, gibbs_burn_in = 1000,
                     is_parallel = TRUE, no_cores = NULL, is_verbose = FALSE){

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
    if (no_cores > K){
      no_cores <- K
    }
    # Create cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)
  }

  # Extract number of observations
  N <- length(x)
  # Number of coefficients
  M <- basis$M + 1

  # Initialize parameters
  if (is.null(dir_a)){
    dir_a <- rep(1 / K, K)
  }
  if (is.null(w_0_mean)){
    w_0_mean <- rep(0, M)
  }
  if (is.null(w_0_cov)){
    w_0_cov <- diag(1, M)
  }
  # Invert the matrix
  prec_0 <- solve(w_0_cov)

  # Matrices for storing results
  weighted_pdf <- matrix(0, nrow = N, ncol = K)  # Store weighted PDFs
  post_prob <- matrix(0, nrow = N, ncol = K)  # Hold responsibilities
  NLL <- vector(mode = "numeric")  # Hold NLL for all MCMC iterations
  C_n <- matrix(0, nrow = N, ncol = K)  # Mixture components
  C_matrix <- matrix(0, nrow = N, ncol = K)  # Total Mixture components

  # Mixing Proportions
  pi_draws <- matrix(0, nrow = gibbs_nsim, ncol = K)
  pi_draws[1, ] <- pi_k

  # Array for storing the coefficient draws for each cluster
  w_draws <- array(data = NA, dim = c(gibbs_nsim, M, K))

  for (k in 1:K){
    w_draws[1, ,k] <- w[, k]
  }

  if (is_parallel){
    # Create design matrix for each observation
    des_mat <- parallel::mclapply(X   = x,
                                  FUN = function(y)
                                    design_matrix(x = basis, obs = y[ ,1])$H,
                                  mc.cores = no_cores)
  }else{
    # Create design matrix for each observation
    des_mat <- lapply(X   = x,
                      FUN = function(y)
                        design_matrix(x = basis, obs = y[ ,1])$H)
  }

  ## ----------------------------------------------------------------------
  # Auxiliary variable model parameters
  ext_des_mat <- list()
  data_y <- list()
  # N1 in 1 column, N0 in second column
  trial_success_mat <- matrix(0, ncol = 2, nrow = N)

  # Iterate over each region
  for (i in 1:N){
    # Total number of reads for each CpG
    N_i <- x[[i]][, 2]
    # Corresponding number of methylated reads for each CpG
    m_i <- x[[i]][, 3]

    # Create extended vector y of length (J x 1)
    y <- vector(mode = "integer")
    for (j in 1:NROW(x[[i]])){
      y <- c(y, rep(1, m_i[j]), rep(0, N_i[j] - m_i[j]))
    }
    data_y[[i]] <- y

    # Col1: Number of successes
    # Col2: Number of failures
    trial_success_mat[i, ] <- c(sum(y), sum(N_i) - sum(y))

    # Create extended design matrix H of dimension (J x M)
    ext_des_mat[[i]] <- as.matrix(des_mat[[i]][rep(1:NROW(des_mat[[i]]),
                                                   N_i), ])
  }

  # TODO: Keep only one design matrix PROBABLE MEMORY ISSUE

  # Run Gibbs sampling
  for (t in 2:gibbs_nsim){

    ## -----------------------------------------------------------------------
    # Compute weighted pdfs for each cluster
    for (k in 1:K){
      # For each element in x, evaluate the BPR log likelihood
      weighted_pdf[ ,k] <- vapply(X   = 1:N,
                                  FUN = function(y)
                                    bpr_likelihood(w = w[, k],
                                                   H = des_mat[[y]],
                                                   data = x[[y]][ ,2:3],
                                                   is_NLL = FALSE),
                                  FUN.VALUE = numeric(1),
                                  USE.NAMES = FALSE)
      weighted_pdf[ ,k] <- log(pi_k[k]) + weighted_pdf[, k]
    }

    # Calculate probabilities using the logSumExp trick for numerical stability
    Z <- apply(weighted_pdf, 1, log_sum_exp)
    # Get actual posterior probabilities, i.e. responsibilities
    post_prob <- exp(weighted_pdf - Z)
    # Evaluate and store the NLL
    NLL  <- c(NLL, (-1) * sum(Z))

    ## ---------------------------------------------------------------------
    # Draw mixture components for ith simulation
    for (i in 1:N){ # Sample one point from a multinomial i.e. ~ Discrete
      C_n[i, ] <- rmultinom(n = 1, size = 1, post_prob[i, ])
    }
    C_matrix <- C_matrix + C_n

    ## ---------------------------------------------------------------------
    # Calculate component counts of each cluster
    N_k <- colSums(C_n)
    # Update mixing proportions using new cluster component counts
    # by sampling from Dirichlet
    pi_k <- as.vector(rdirichlet(n = 1, alpha = dir_a + N_k))
    pi_draws[t, ] <- pi_k


    ## --------------------------------------------------------------------
    for (k in 1:K){
      # Which regions are assigned to cluster k
      C_k_idx <- which(C_n[, k] == 1)

      # Concatenate data from all regions in cluster k
      H <- do.call(rbind, ext_des_mat[C_k_idx])

      # Concatenate y from all regions in cluster k
      y <- do.call(c, data_y[C_k_idx])

      # Add all successes and failures from all regions in cluster k
      N1_N0 <- colSums(trial_success_mat[C_k_idx, ])

      # Compute posterior variance of w
      V <- solve(prec_0 + crossprod(H, H))

      # Update Mean of z
      mu_z <- H %*% w[, k]
      # Draw latent variable z from its full conditional: z | \w, y, X
      z <- rep(0, sum(N1_N0))
      z[y == 1] <- rtruncnorm(N1_N0[1], mean = mu_z[y == 1], sd = 1,
                              a = 0, b = Inf)
      z[y == 0] <- rtruncnorm(N1_N0[2], mean = mu_z[y == 0], sd = 1,
                              a = -Inf, b = 0)

      # Compute posterior mean of w
      Mu <- V %*% (prec_0 %*% w_0_mean + crossprod(H, z))
      # Draw variable \w from its full conditional: \w | z, X
      w[, k] <- c(rmvnorm(1, Mu, V))

      # Store the \theta draws
      w_draws[t, , k] <- w[, k]
    }


#     # Update basis function coefficient vector w for each cluster
#     # If parallel mode is ON
#     if (is_parallel){
#       # Parallel optimization for each cluster k
#       w <- foreach::"%dopar%"(obj = foreach::foreach(k = 1:K,
#                                                      .combine = cbind),
#                               ex  = {
#                                 out <- gibbs_bpr(par       = w[ ,k],
#                                              fn        = sum_weighted_bpr_lik,
#                                              gr        = sum_weighted_bpr_grad,
#                                              method    = opt_method,
#                                              control   = list(maxit = opt_itnmax),
#                                              x         = x,
#                                              des_mat   = des_mat,
#                                              post_prob = post_prob[ ,k],
#                                              is_NLL    = TRUE)
#                               })
#     }else{
#       # Sequential optimization for each clustrer k
#       w <- foreach::"%do%"(obj = foreach::foreach(k = 1:K,
#                                                   .combine = cbind),
#                            ex  = {
#                              out <- optim(par       = w[ ,k],
#                                           fn        = sum_weighted_bpr_lik,
#                                           gr        = sum_weighted_bpr_grad,
#                                           method    = opt_method,
#                                           control   = list(maxit = opt_itnmax),
#                                           x         = x,
#                                           des_mat   = des_mat,
#                                           post_prob = post_prob[ ,k],
#                                           is_NLL    = TRUE)$par
#                            })
#     }
  }

  if (is_parallel){
    # Stop parallel execution
    parallel::stopCluster(cl)
  }

  obj <- structure(list(K = K,
                        N = N,
                        w = w,
                        pi_k = pi_k,
                        NLL = NLL,
                        basis = basis,
                        C_matrix = C_matrix,
                        w_draws = w_draws,
                        pi_draws = pi_draws,
                        post_prob = post_prob),
                   class = "bpr_fdmm")
  return(obj)
}
