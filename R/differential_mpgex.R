#' Predict differntial gene expression from differential methylation
#'
#' \code{differential_mpgex} is a function that wraps all the necessary
#' subroutines for performing predictions on differential gene expression.
#' Initially, it optimizes the parameters of the basis functions so as to
#' learn the methylation profiles for the control case and the treatment case.
#' Then the two learned profiles for each promoter region are subtracted and
#' new coefficients showing the difference between the two profiles are learnt
#' using the Basis Linear Model. These coefficients are given as input features
#' for performing linear regression in order to predict/regress the
#' corresponding differential gene expression data.
#'
#' @inheritParams wrapper_mpgex
#' @param X The binomial distributed observations. A list containing two lists
#'  for control and treatment samples. Each of the two lists has a nested list
#'  where each element is an L x 3 dimensional matrix.
#' @param lambda Regularization term when performing Basis Linear Model fitting
#'
#' @return An mpgex object consisting of the following elements: ....
#'
#' @seealso \code{\link{bpr_optim}}, \code{\link{bpr_likelihood}},
#'  \code{\link{blm}}, \code{\link{design_matrix}}
#'
#' @examples
#' obs <- list(control = bpr_control_data, treatment = bpr_treatment_data)
#' Y   <- list(control = gex_control_data, treatment = gex_treatment_data)
#' basis <- rbf.object(M = 5, gamma = 0.3)
#' out   <- differential_mpgex(X = obs, Y = Y, basis = basis, lambda = 1e-02)
#'
#' @export
differential_mpgex <- function(formula = NULL, X, Y, train_ind = NULL,
                         basis = NULL, w = NULL, train_perc = 0.7,
                         method = "CG", itnmax = 100, lambda = 0){

  # Optimize the BPR function for control samples
  out_contr_opt <- bpr_optim(x      = X$control,
                             w      = w,
                             basis  = basis,
                             method = method,
                             itnmax = itnmax)

  # Optimize the BPR function for treatment samples
  out_treat_opt <- bpr_optim(x       = X$treatment,
                             w       = out_contr_opt$w,
                             basis   = out_contr_opt$basis,
                             method  = method,
                             itnmax  = itnmax)


  out_diff_meth <- fit_diff_meth(control    = out_contr_opt,
                                 treatment  = out_treat_opt,
                                 diff_basis = out_contr_opt$basis,
                                 lambda     = lambda)

  diff_expr <- (Y$control + 1e-01) / (Y$treatment + 1e-01)

  # Create training and test sets
  dataset <- partition_data(X = out_diff_meth$W_opt,
                            Y = diff_expr,
                            train_ind  = train_ind,
                            train_perc = train_perc)

  # Predict gene expression from methylation profiles
  mpgex <- predict_gex(formula = formula,
                       train   = dataset$train,
                       test    = dataset$test)

  # Create a 'diff_mpgex' object
  mpgex$basis       <- out_contr_opt$basis
  mpgex$diff_basis  <- out_diff_meth$diff_basis
  mpgex$W_opt       <- out_diff_meth$W_opt
  mpgex$diff_meth   <- out_diff_meth$diff_meth
  mpgex$method <- method
  mpgex$itnmax <- itnmax
  mpgex$train  <- dataset$train
  mpgex$test   <- dataset$test
  mpgex$contr_opt <- out_contr_opt
  mpgex$treat_opt <- out_treat_opt
  mpgex$diff_expr <- diff_expr
  mpgex$train_ind <- train_ind

  class(mpgex) <- "diff_mpgex"

  return(mpgex)
}
