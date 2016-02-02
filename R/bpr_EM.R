#' EM algorithm for BPR mixture model
#'
#' \code{bpr_EM} implements the EM algorithm for performing clustering on DNA
#'  methylation profiles, where the observation model is the Binomial
#'  distributed Probit Regression function, \code{\link{bpr_likelihood}}.
#'
#' @param x A list of elements of length N, where each element is an L x 3
#' matrix of observations, where 1st column contains the locations. The 2nd
#' and 3rd columns contain the total trials and number of successes at the
#' corresponding locations, repsectively.
#' @param K Integer denoting the number of clusters K.
#' @param pi_k Vector of length K, denoting the mixing proportions.
#' @param w A MxK matrix, where each column contains the basis function
#' coefficients for the corresponding cluster.
#' @param em_max_iter Integer denoting the maximum number of EM iterations.
#' @param epsilon_conv Numeric denoting the convergence parameter for EM.
#' @param is_verbose Logical, print results during EM iterations
#' @inheritParams bpr_optim.matrix
#'
#' @examples
#' ex_data <- bpr_data
#' out_opt <- bpr_EM(x = ex_data, em_max_iter = 10, opt_itnmax = 10)
#'
#' @export
bpr_EM <- function(x, K = 2, pi_k = NULL, w = NULL, basis = NULL,
                   em_max_iter = 100, epsilon_conv = 1e-05, opt_method = "CG",
                                        opt_itnmax = 100, is_verbose = FALSE){
  # Check that x is a list object
  assertthat::assert_that(is.list(x))

  # Extract number of observations
  N <- length(x)
  assertthat::assert_that(N > 0)

  # Store weighted PDFs
  weighted_pdf <- matrix(0, nrow = N, ncol = K)

  # Initialize and store NLL for each EM iteration
  NLL <- c(1e+40)

  # Perform checks for initial parameter values
  out <- do_EM_checks(x = x,
                      K = K,
                      pi_k = pi_k,
                      w = w,
                      basis = basis,
                      opt_method = opt_method,
                      opt_itnmax = opt_itnmax)
  w   <- out$w
  basis <- out$basis
  pi_k <- out$pi_k

  # Create design matrix for each observation
  des_mat <- lapply(X   = x,
                    FUN = function(y)
                      design_matrix(x = basis, obs = y[ ,1]))

  # Run EM algorithm until convergence
  for (t in 1:em_max_iter){

    #
    # E-Step -----------------------------------------------
    #
    # Compute weighted pdfs for each cluster
    for (k in 1:K){
      # For each element in x, evaluate the BPR log likelihood
      weighted_pdf[ ,k] <- vapply(X   = 1:N,
                                  FUN = function(y)
                                    bpr_likelihood(w = w[ ,k],
                                                   H = des_mat[[y]]$H,
                                                   data = x[[y]][ ,2:3],
                                                   is_NLL = FALSE),
                                  FUN.VALUE = numeric(1),
                                  USE.NAMES = FALSE)
      weighted_pdf[ ,k] <- log(pi_k[k]) + weighted_pdf[ ,k]
    }
    # Calculate probabilities using the logSumExp trick for numerical stability
    Z <- apply(weighted_pdf, 1, log_sum_exp)

    # Get actual posterior probabilities, i.e. responsibilities
    post_prob <- exp(weighted_pdf - Z)

    # Evaluate and store the NLL
    NLL  <- c(NLL, (-1) * sum(Z))

    #
    # M-Step -----------------------------------------------
    #
    # Compute sum of posterior probabilities for each cluster
    N_k <- colSums(post_prob)

    # Update mixing proportions for each cluster
    pi_k <- N_k / N

    # Update basis function coefficient vector w for each cluster
    for (k in 1:K){
      w[ ,k] <- optim(par       = w[ ,k],
                      fn        = sum_weighted_bpr_lik,
                      gr        = sum_weighted_bpr_grad,
                      method    = opt_method,
                      control   = list(maxit = opt_itnmax),
                      x         = x,
                      des_mat   = des_mat,
                      post_prob = post_prob[ ,k],
                      is_NLL    = TRUE)$par
    }

    if (is_verbose){
      cat("It:\t",t, "\tNLL:\t", NLL[t+1],
              "\tNLL_diff:\t", NLL[t] - NLL[t + 1], "\n")
    }

    if (NLL[t + 1] > NLL[t]){
      stop("Negative Log Likelihood increases - Stopping EM!\n")
    }

    # Check for convergence
    if (NLL[t] - NLL[t + 1] < epsilon_conv){
      break
    }
  }

  # Check if EM converged in the given maximum iterations
  if (t == em_max_iter){
    warning("EM did not converge with the given maximum iterations!\n")
  }

  # Add names to the estimated parameters for clarity
  names(pi_k) <- paste0("clust", 1:K)
  colnames(w) <- paste0("clust", 1:K)

  # Get hard cluster assignments for each observation
  clust_assign <- apply(X = post_prob,
                        MARGIN = 1,
                        FUN = function(x) which(x == max(x, na.rm = TRUE)))

  # Perform model selection
  total_params <- (K - 1) + K * NROW(w)

  # Bayesian Information Criterion
  BIC <- 2 * NLL + total_params * log(N)

  # Akaike Iformation Criterion
  AIC <- 2 * NLL + 2 * total_params

  # Integrated Complete Likelihood criterion
  entropy <- (-1) * sum(post_prob * log(post_prob), na.rm = TRUE)
  ICL <- BIC + entropy

  obj <- structure(list(K = K,
                        N = N,
                        w = w,
                        pi_k = pi_k,
                        basis = basis,
                        post_prob = post_prob,
                        clust_assign = clust_assign,
                        NLL = NLL,
                        BIC = BIC,
                        AIC = AIC,
                        ICL = ICL),
                   class = "mpgex_EM")
  return(obj)
}


# Internal function to make all the appropriate type checks.
do_EM_checks <- function(x, K = 2, pi_k = NULL,  w = NULL, basis = NULL,
                                        opt_method = "CG", opt_itnmax = 100){
  if (is.null(basis)){
    basis <- polynomial.object()
  }
  if (is.null(w)){
    w <- rep(0.2, basis$M + 1)

    # Optimize the BPR function for each element in x
    out_opt <- bpr_optim(x           = x,
                         w           = w,
                         basis       = basis,
                         method      = opt_method,
                         itnmax      = opt_itnmax)

    # Keep only the optimized coefficients
    W_opt <- out_opt$W_opt

    # Use Kmeans with random starts
    cl <- kmeans(W_opt, K, nstart = 25)
    # Get the mixture components
    C_n <- cl$cluster
    # Mean for each cluster
    w <- t(cl$centers)

    # Mixing proportions
    if (is.null(pi_k)){
      N <- length(x)
      pi_k <- as.vector(table(C_n) / N)
    }
  }
  if (NROW(w) != (basis$M + 1) ){
    stop("Coefficients vector should be M+1, M: number of basis functions!")
  }
  return(list(w = w, basis = basis, pi_k = pi_k))
}
