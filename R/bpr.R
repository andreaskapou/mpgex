#' BPR log likelihood function
#'
#' \code{bpr_likelihood} evaluates the Binomial distributed Probit regression
#'  log likelihood function for a given set of coefficients, observations and
#'  a design matrix.
#'
#' @section Mathematical formula:
#'  The Binomial distributed Probit regression log likelihood function
#'  is computed by the following formula:
#'  \deqn{log p(y | f, w) = \sum_{l=1}^{L} log Binom(m_{l} | t_{l}, \Phi(w^{t}h(x_{l})))}
#'  where h(x_{l}) are the basis functions.
#'
#' @param w A vector of parameters (i.e. coefficients of the basis functions)
#' @param H The \code{L x M} matrix design matrix, where L is the number
#'  of observations and M the number of basis functions.
#' @param data An \code{L x 2} matrix containing in the 1st column the total
#'  number of trials and in the 2nd the number of successes. Each row
#'  corresponds to each row of the design matrix.
#' @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#'  returned.
#'
#' @return the log likelihood
#'
#' @seealso \code{\link{bpr_gradient}}, \code{\link{design_matrix}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' obs <- c(0,.2,.5, 0.6)
#' des_mat <- design_matrix(obj, obs)
#' H <- des_mat$H
#' w <- c(.1,.1,.1)
#' data <- matrix(c(10,12,15,7,9,8), ncol=2)
#' lik <- bpr_likelihood(w, H, data)
#'
#' @importFrom stats pnorm dbinom
#' @export
bpr_likelihood <- function(w, H, data, is_NLL = FALSE){
  total <- data[ ,1]
  succ  <- data[ ,2]

  # Predictions of the target variables
  # Compute the cdf of N(0,1) distribution (i.e. probit function)
  Phi <- pnorm(H %*% w)

  # In extreme cases where probit is 0 or 1, subtract a tiny number
  # so we can evaluate the log(0) when computing the Binomial
  Phi[which(Phi > (1 - 1e-289))] <- 1 - 1e-289
  Phi[which(Phi < 1e-289)] <- 1e-289

  # Compute the log likelihood
  res <- sum(dbinom(x = succ, size = total, prob = Phi, log = TRUE)) -
                                                       1/2 * t(w) %*% w

  # If we required the Negative Log Likelihood
  if (is_NLL){
    res <- (-1) * res
  }
  return(res)
}


#' Gradient of the BPR log likelihood function
#'
#' \code{bpr_gradient} computes the gradient w.r.t the coefficients w of
#'  the Binomial distributed Probit regression log likelihood function.
#'
#' @section Mathematical formula:
#'  The gradient of the Binomial distributed Probit regression log likelihood
#'  function w.r.t to w is computed by the following formula:
#' \deqn{log p(y | f, w) = \sum_{l=1}^{L} log Binom(m_{l} | t_{l}, \Phi(w^{t}h(x_{l})))}
#'
#' @inheritParams bpr_likelihood
#'
#' @return the gradient vector of the log likelihood w.r.t the vector of
#' coefficients w
#'
#' @seealso \code{\link{bpr_likelihood}}, \code{\link{design_matrix}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' obs <- c(0,.2,.5)
#' des_mat <- design_matrix(obj, obs)
#' H <- des_mat$H
#' w <- c(.1,.1,.1)
#' data <- matrix(c(10,12,15,7,9,8), ncol=2)
#' gr <- bpr_gradient(w, H, data)
#'
#' @importFrom stats pnorm dnorm
#' @export
bpr_gradient <- function(w, H, data, is_NLL = FALSE){
  total <- data[ ,1]
  succ  <- data[ ,2]

  # Predictions of the target variables
  g <- as.vector(H %*% w)
  # Compute the cdf of N(0,1) distribution (i.e. probit function)
  Phi <- pnorm(g)

  # In extreme cases where probit is 0 or 1, subtract a tiny number
  # so we can evaluate the log(0) when computing the Binomial
  Phi[which(Phi > (1 - 1e-289))] <- 1 - 1e-289
  Phi[which(Phi < 1e-289)] <- 1e-289

  # Compute the density of a N(0,1) distribution
  N <- dnorm(g)
  N[which(N < 1e-289)] <- 1e-289

  # Compute the gradient vector w.r.t the coefficients w
  gr <- (N * (succ * (1/Phi) - (total - succ) * (1 / (1 - Phi)))) %*% H - w

#   if (NROW(H) == 1){
#     gr <- (succ * (1 / Phi) - (total - succ) * (1 / (1 - Phi))) * N %*% H - w
#   }else{
#     gr <- (t(succ) %*% diag(1 / Phi) - t(total - succ) %*%
#              diag(1 / (1 - Phi))) %*% diag(N) %*% H - w
#
#     (N * (succ * (1/Phi) - (total - succ) * (1 / (1 - Phi)))) %*% H - w
#   }

  # If we required the Negative Log Likelihood
  if (is_NLL){
    gr <- (-1) * gr
  }
  return(gr)
}


#' Sum of weighted BPR log likelihoods
#'
#' \code{sum_weighted_bpr_lik} computes the sum of the BPR log likelihoods for
#' each elements of x, and then weights them by the corresponding posterior
#' probabilities. This function is mainly used for the M-Step of the EM
#' algorithm \code{\link{bpr_EM}}.
#'
#' @param w A vector of parameters (i.e. coefficients of the basis functions)
#' @param x A list of elements of length N, where each element is an L x 3
#'  matrix of observations, where 1st column contains the locations. The 2nd
#'  and 3rd columns contain the total trials and number of successes at the
#'  corresponding locations, repsectively.
#' @param des_mat A list of length N, where each element contains the
#'  \code{L x M} design matrices, where L is the number of observations and M
#'  the number of basis functions.
#' @param post_prob A vector of length N containing the posterior probabilities
#'  fore each element of list x, respectively.
#' @param is_NLL Logical, indicating if the Negative Log Likelihood should be
#'  returned.
#'
#' @return The weighted sum of BPR log likelihoods
#'
#' @seealso \code{\link{bpr_likelihood}}, \code{\link{bpr_EM}}
#'
#' @export
sum_weighted_bpr_lik <- function(w, x, des_mat, post_prob, is_NLL = TRUE){
  N <- length(x)

  # TODO: Create tests
  # For each element in x, evaluate the BPR log likelihood
  res <- vapply(X   = 1:N,
                FUN = function(y) bpr_likelihood(w = w,
                                                 H = des_mat[[y]]$H,
                                                 data = x[[y]][ ,2:3],
                                                 is_NLL = is_NLL),
                FUN.VALUE = numeric(1),
                USE.NAMES = FALSE)

  # Return the dot product of the result and the posterior probabilities
  return(post_prob %*% res)
}


#' Sum of weighted gradients of the BPR log likelihood
#'
#' \code{sum_weighted_bpr_grad} computes the sum of the gradients of BPR log
#' likelihood for each elements of x, and then weights them by the
#' corresponding posterior probabilities. This function is mainly used for the
#' M-Step of the EM algorithm \code{\link{bpr_EM}}.
#'
#' @inheritParams sum_weighted_bpr_lik
#'
#' @return A vector with weighted sum of the gradients of BPR log likelihood.
#'
#' @seealso \code{\link{bpr_gradient}}, \code{\link{bpr_EM}}
#'
#' @export
sum_weighted_bpr_grad <- function(w, x, des_mat, post_prob, is_NLL = TRUE){
  N <- length(x)

  # TODO: Create tests
  # For each element in x, evaluate the gradient of the BPR log likelihood
  res <- vapply(X   = 1:N,
                FUN = function(y) bpr_gradient(w = w,
                                               H = des_mat[[y]]$H,
                                               data = x[[y]][ ,2:3],
                                               is_NLL = is_NLL),
                FUN.VALUE = numeric(length(w)),
                USE.NAMES = FALSE)

  # Return the dot product of the result and the posterior probabilities
  return(post_prob %*% t(res))
}
