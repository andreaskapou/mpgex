#' Predict differntial gene expression from differential methylation
#'
#' \code{mpgex_diff_regr} is a function that wraps all the necessary
#' subroutines for performing predictions on differential gene expression.
#' Initially, it optimizes the parameters of the basis functions so as to
#' learn the methylation profiles for the control case and the treatment case.
#' Then the two learned profiles for each promoter region are subtracted and
#' new coefficients showing the difference between the two profiles are learned
#' using the Basis Linear Model. These coefficients are given as input features
#' for performing linear regression in order to predict/regress the
#' corresponding differential gene expression data.
#'
#' @param x The binomial distributed observations. A list containing two lists
#'  for control and treatment samples. Each of the two lists has a nested list
#'  where each element is an L x 3 dimensional matrix.
#' @param y The gene expression data. A list containing two vectors for control
#' and treatment samples.
#' @param diff_basis The differential basis object, if NULL, the main basis
#'  object will be created.
#' @param lambda Regularization term when performing Basis Linear Model fitting
#' @param x_eval Integer denoting the number of evaluation points when
#'  learning the differential methylation profiles.
#' @inheritParams mpgex_regr
#'
#' @return An mpgex object consisting of the following elements:
#'
#' @seealso \code{\link{bpr_optim}}, \code{\link{bpr_likelihood}},
#'  \code{\link{blm}}, \code{\link{polynomial.object}},
#'  \code{\link{rbf.object}}
#'
#' @examples
#' obs <- list(control = bpr_control_data, treatment = bpr_treatment_data)
#' y   <- list(control = gex_control_data, treatment = gex_treatment_data)
#' basis <- rbf.object(M = 3)
#' set.seed(1234)
#' out   <- mpgex_diff_regr(x = obs, y = y, basis = basis, is_parallel = FALSE,
#'                          opt_itnmax = 50, lambda = 1e-02)
#'
#' @export
mpgex_diff_regr <- function(formula = NULL, x, y, model_name = "svm", w = NULL,
                            basis = NULL, diff_basis = NULL, train_ind = NULL,
                            train_perc = 0.7,
                            fit_feature = NULL, opt_method = "CG",
                            opt_itnmax = 500, lambda = 0, is_parallel = TRUE,
                            no_cores = NULL, x_eval = 50, is_summary = TRUE){

  # Learn methylation profiles for control samples
  message("Learning control methylation profiles ...\n")
  out_contr_opt <- bpr_optim(x           = x$control,
                             w           = w,
                             basis       = basis,
                             fit_feature = NULL,
                             opt_method  = opt_method,
                             opt_itnmax  = opt_itnmax,
                             is_parallel = is_parallel,
                             no_cores    = no_cores)

  # Learn methylation profiles for treatment samples
  message("Learning treatment methylation profiles ...\n")
  out_treat_opt <- bpr_optim(x           = x$treatment,
                             w           = out_contr_opt$w,
                             basis       = out_contr_opt$basis,
                             fit_feature = NULL,
                             opt_method  = opt_method,
                             opt_itnmax  = opt_itnmax,
                             is_parallel = is_parallel,
                             no_cores    = no_cores)

  # Choose by default an RBF basis with 3 kernels
  if (is.null(diff_basis)){
    diff_basis <- out_contr_opt$basis
  }

  # Learn differential methylation profile from control and treatment samples
  message("Learning differential methylation profiles ...\n")
  out_diff_meth <- learn_diff_meth(control     = out_contr_opt,
                                   treatment   = out_treat_opt,
                                   diff_basis  = diff_basis,
                                   fit_feature = fit_feature,
                                   lambda      = lambda,
                                   x_eval      = x_eval,
                                   is_parallel = is_parallel,
                                   no_cores    = no_cores)

  # TODO: Should this be division or subtraction or something more complicated
  # Compute differential gene expression
  diff_expr <- (y$control + 1e-02) - (y$treatment + 1e-02)

  # Create training and test sets
  message("Partitioning to test and train data ...\n")
  dataset <- partition_data(x          = out_diff_meth$W_opt,
                            y          = diff_expr,
                            train_ind  = train_ind,
                            train_perc = train_perc)

  # Train regression model from methylation profiles
  message("Training linear regression model ...\n")
  train_model <- train_model_gex(formula    = formula,
                                 model_name = model_name,
                                 train      = dataset$train,
                                 is_summary = is_summary)

  # Predict gene expression from methylation profiles
  message("Making predictions ...\n")
  predictions <- predict_model_gex(model      = train_model$gex_model,
                                   test       = dataset$test,
                                   is_summary = is_summary)
  message("Done!\n\n")

  # Create 'mpgex_diff_regr' object
  obj <- structure(list(formula      = formula,
                        model_name   = model_name,
                        opt_method   = opt_method,
                        opt_itnmax   = opt_itnmax,
                        train_ind    = dataset$train_ind,
                        gex_model    = train_model$gex_model,
                        train_pred   = train_model$train_pred,
                        test_pred    = predictions$test_pred,
                        train_errors = train_model$train_errors,
                        test_errors  = predictions$test_errors,
                        diff_expr    = diff_expr,
                        diff_meth    = out_diff_meth$diff_meth,
                        diff_basis   = diff_basis,
                        xs_mat       = out_diff_meth$xs_mat,
                        train        = dataset$train,
                        test         = dataset$test,
                        basis        = out_contr_opt$basis,
                        W_opt        = out_diff_meth$W_opt,
                        Mus          = out_diff_meth$Mus),
                   class = "mpgex_diff_regr")
  return(obj)
}
