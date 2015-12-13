#' Predict differntial gene expression from scaled differential methylation
#'
#' \code{scaled_mpgex} is a function that wraps all the necessary
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
#'  regression. If NULL, the simple linear regression model is used.
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
#' obs <- list(control = bpr_control_data, treatment = bpr_treatment_data)
#' Y   <- list(control = gex_control_data, treatment = gex_treatment_data)
#' basis <- rbf.object(M = 5, gamma = 0.3)
#' out   <- scaled_mpgex(X = obs, Y = Y, basis = basis)
#'
#' @export
scaled_mpgex <- function(formula = NULL, X, Y, train_ind = NULL,
                               basis = NULL, w = NULL, train_perc = 0.7,
                                            method = "CG", itnmax = 100){

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

  W_opt <- NULL
  W_opt <- out_treat_opt$W_opt / out_contr_opt$W_opt
  diff_expr <- Y$control - Y$treatment

  # Create training and test sets
  dataset <- partition_data(X = W_opt,
                            Y = diff_expr,
                            train_ind  = train_ind,
                            train_perc = train_perc)

  # Predict gene expression from methylation profiles
  mpgex <- predict_gex(formula = formula,
                       train   = dataset$train,
                       test    = dataset$test)

  # Create an 'mpgex' object
  mpgex$basis  <- out_treat_opt$basis
  mpgex$W_opt  <- W_opt
  mpgex$method <- method
  mpgex$itnmax <- itnmax
  mpgex$contr_opt <- out_contr_opt
  mpgex$treat_opt <- out_treat_opt
  mpgex$diff_expr <- diff_expr
  mpgex$train_ind <- train_ind

  class(mpgex) <- "diff-mpgex"

  return(mpgex)
}
