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
#' @param data An \code{L x 2} matrix containing in the 1st column the number of
#'  total reads and and in the 2nd the number of methylated reads. Each row
#'  corresponds to each row of the design matrix.
#'
#' @return the log likelihood
#'
#' @seealso \code{\link{bpr_derivative}}, \code{\link{design_matrix}}
#'
#' @examples
#' obj <- polynomial.object(M=2)
#' obs <- c(0,.2,.5)
#' des_mat <- design_matrix(obj, obs)
#' H <- des_mat$H
#' w <- c(.1,.1,.1)
#' data <- matrix(c(10,12,15,7,9,8), ncol=2)
#' lik <- bpr_likelihood(w, H, data)
#'
#' @export
bpr_likelihood <- function(w, H, data){
  total <- data[ ,1]
  meth  <- data[ ,2]

  # Predictions of the target variables
  g <- as.vector(H %*% w)
  # Compute the cdf of N(0,1) distribution (i.e. probit function)
  Phi <- as.vector(pnorm(g))

  # Compute the log likelihood
  res <- sum(dbinom(x = meth, size = total, prob = Phi, log = TRUE))
  return(res)
}


#' Derivative of the BPR log likelihood function
#'
#' \code{bpr_derivative} computes the derivative w.r.t the coefficients w of
#'  the Binomial distributed Probit regression log likelihood function.
#'
#' @section Mathematical formula:
#'  The derivative of the Binomial distributed Probit regression log likelihood
#'  function w.r.t to w is computed by the following formula:
#' \deqn{log p(y | f, w) = \sum_{l=1}^{L} log Binom(m_{l} | t_{l}, \Phi(w^{t}h(x_{l})))}
#'
#' @inheritParams bpr_likelihood
#'
#' @return the derivative vector of the log likelihood w.r.t the vector of
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
#' der <- bpr_derivative(w, H, data)
#'
#' @export
bpr_derivative <- function(w, H, data){
  total <- data[ ,1]
  meth  <- data[ ,2]

  # Predictions of the target variables
  g <- as.vector(H %*% w)
  # Compute the cdf of N(0,1) distribution (i.e. probit function)
  Phi <- as.vector(pnorm(g))
  # Compute the density of a N(0,1) distribution
  N <- as.vector(dnorm(g))

  # Compute the derivative vector w.r.t the coefficients w
  der <- (t(meth) %*% diag(1 / Phi) - t(total - meth) %*%
            diag(1 / (1 - Phi))) %*% diag(N) %*% H

  return(der)
}
