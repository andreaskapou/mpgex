#' Predict gene expression from methylation profiles
#'
#' \code{wrapper_mpgex} is a function that wraps all the necessary subroutines
#' for performing predictions on gene expressions. Initially, it optimizes the
#' parameters of the basis functions so as to learn the methylation profiles.
#' Then uses the learned parameters / coefficients of the basis functions as
#' input features for performing linear regression in order to predict/regress
#' the corresponding gene expression data.
#'
#' @param formula An object of class \code{\link[stats]{formula}} needed when
#'  calling the \code{\link[stats]{lm}} function for performing linear
#'  regression. If NULL, the simple linear regression method is used.
#' @param X The binomial distributed observations, which has to be a list
#'  where each element is an L x 3 dimensional matrix.
#' @param Y Corresponding gene expression data for each element of the list X
#' @param train_ind Optional argument containing the indices for splitting the
#'  dataset to train and test sets.
#' @param basis Optional basis function object, default is
#'  \code{\link{polynomial.object}}
#' @param w Optional vector of initial parameter / coefficient values.
#' @param train_perc Optional parameter for defining the percentage of the
#'  dataset to be used for training set, the remaining will be the test set.
#' @param method Parameter for defining the method to be used in the
#'  optimization procedure, see \code{\link[stats]{optim}}.
#' @param itnmax Optional parameter for defining the max number of iterations
#'  of the optimization procedure, see \code{\link[stats]{optim}}.
#'
#' @return An mpgex object consisting of the following elements: ...
#'
#' @seealso \code{\link{bpr_optim}}, \code{\link{bpr_likelihood}},
#'  \code{\link{polynomial.object}}, \code{\link{rbf.object}},
#'  \code{\link{design_matrix}}
#'
#' @examples
#' obs <- bpr_data
#' Y   <- gex_data
#' basis <- rbf.object(M = 5, gamma = 0.3)
#' out   <- wrapper_mpgex(X = obs, Y = Y, basis = basis)
#'
#' @export
wrapper_mpgex <- function(formula = NULL, X, Y, model_name = "lm",
                          train_ind = NULL, basis = NULL, w = NULL,
                          train_perc = 0.7, method = "CG", itnmax = 100){

  # Optimize the BPR function for each element in x
  out_opt <- bpr_optim(x = X,
                       w = w,
                       basis = basis,
                       method = method,
                       itnmax = itnmax)

  # Create training and test sets
  dataset <- partition_data(X = out_opt$W_opt,
                            Y = Y,
                            train_ind  = train_ind,
                            train_perc = train_perc)

  # Predict gene expression from methylation profiles
  mpgex <- predict_gex(formula    = formula,
                       model_name = model_name,
                       train      = dataset$train,
                       test       = dataset$test)

  # Create an 'mpgex' object
  mpgex$basis  <- out_opt$basis
  mpgex$W_opt  <- out_opt$W_opt
  mpgex$Mus    <- out_opt$Mus
  mpgex$train  <- dataset$train
  mpgex$test   <- dataset$test
  mpgex$method <- method
  mpgex$itnmax <- itnmax
  mpgex$train_ind <- dataset$train_ind
  class(mpgex) <- "mpgex"

  return(mpgex)
}
