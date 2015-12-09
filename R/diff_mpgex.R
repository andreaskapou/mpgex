#' Predict differntial gene expression from differential methylation
#'
#' \code{differential_mpgex} is a function that wraps all the necessary
#' subroutines for performing predictions on differential gene expression.
#' Initially, it optimizes the parameters of the basis functions so as to
#' learn the methylation profiles for the control case. Then uses these
#' parameters as constants and learns the methylation profiles of the treatment
#' as a scaled version of the control profiles. The coefficients showing
#' the difference between the two profiles is given as input features for
#' performing linear regression in order to predict/regress
#' the corresponding differential gene expression data.
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
#' @return An mpgex object consisting of the following elements:
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
differential_mpgex <- function(formula = NULL, X, Y, train_ind = NULL, basis = NULL,
                          w = NULL, train_perc = 0.7, method = "CG", itnmax = 100){

  # -----------------------------------
  # Compute the optimized parameters of the BPR function for control samples
  out_contr_opt <- bpr_optim(x = X$control,
                             w = w,
                             basis = basis,
                             method = method,
                             itnmax = itnmax)

  # -----------------------------------
  # Compute the optimized parameters of the BPR function for treatment samples
  out_treat_opt <- bpr_scaled_optim(x = X$treatment,
                                    W_contr = out_contr_opt$W_opt,
                                    w       = out_contr_opt$w,
                                    basis   = out_contr_opt$basis,
                                    method  = method,
                                    itnmax  = itnmax)

  X_prof <- NULL
  X_prof <- out_treat_opt$W_opt / out_contr_opt$W_opt

  diff_expr <- Y$control - Y$treatment
  # ------------------------------------
  # Create training and test sets
  dataset <- partition_data(X = X_prof,
                            Y = diff_expr,
                            train_ind  = train_ind,
                            train_perc = train_perc)
  train   <- dataset$train
  test    <- dataset$test

  # ------------------------------------
  # Make predictions based on the methylation profiles singatures.
  if (is.null(formula)){
    fit_model  <- stats::lm(formula = Y ~ ., data = train)
  }else{
    fit_model  <- stats::lm(formula = formula, data = train)
  }
  train_pred <- stats::predict(object = fit_model)

  regressors <- 1:(NCOL(test) - 1)
  if (! NCOL(test) == 2){
    test_pred  <- stats::predict(object  = fit_model,
                                 newdata = test[ ,regressors])
  }else{
    test_pred  <- stats::predict(object  = fit_model,
                                 newdata = data.frame(X = test[ ,regressors]))
  }

  # Compute the prediction errors
  train_rmserr <- pracma::rmserr(train$Y, train_pred, summary = TRUE)
  test_rmserr  <- pracma::rmserr(test$Y,  test_pred,  summary = TRUE)

  # Create an 'mpgex' object
  mpgex <- list(X = X,
                control_opt   = out_contr_opt,
                treatment_opt = out_treat_opt,
                diff_expr     = diff_expr,
                X_prof  = X_prof,
                basis   = out_treat_opt$basis,
                method  = method,
                itnmax  = itnmax,
                train_data   = train,
                test_data    = test,
                fit_model    = fit_model,
                train_pred   = train_pred,
                test_pred    = test_pred,
                train_rmserr = train_rmserr,
                test_rmserr  = test_rmserr)
  class(mpgex) <- "diff-mpgex"

  return(mpgex)
}
