#' Learn differential methylation profiles using BLM
#'
#' \code{learn_diff_meth} learns differential methylation profiles and returns
#'  the fitted coefficients using the Basis Linear Model.
#'
#' @param control The optimized values for the control samples.
#' @param treatment The optimized values for the treatment samples.
#' @param diff_basis The basis object for differential methylation.
#' @param fit_feature Additional feature on how well the profile fits the
#'  methylation data.
#' @param lambda The regularization parameter when fitting the BLM.
#' @param x_evals Integer denoting the number of evaluation points.
#'
#' @return The differential methylation profiles and the corresponding fitted
#'  coefficients from the BLM model.
#'
#' @export
learn_diff_meth <- function(control, treatment, diff_basis, fit_feature = NULL,
                                                    lambda = 0, x_evals = 50){
  # Create matrices from lists for efficiency
  contr_W     <- control$W_opt
  treat_W     <- treatment$W_opt
  contr_extr  <- control$x_extrema
  treat_extr  <- treatment$x_extrema

  contr_basis <- control$basis
  treat_basis <- treatment$basis

  # If we have 'rbf' basis, extract the Mus parameters
  if (is(contr_basis, "rbf")){
    if (is.null(contr_basis$mus)){
      contr_Mus <- control$Mus
      treat_Mus <- treatment$Mus
    }
  }

  # Data frame for storing all the coefficients for each element of list x
  if (is.null(fit_feature)){
    num_features <- diff_basis$M + 1
  }else{
    num_features <- diff_basis$M + 2
  }

  # Matrix for storing fitted coefficients
  N <- NROW(contr_W)
  W_opt <- matrix(NA_real_, nrow = N, ncol = num_features)
  colnames(W_opt) <- paste("w", seq(1, num_features), sep = "")

  # Matrix for storing the centers of RBFs if object class is 'rbf'
  Mus <- NULL
  if (is(diff_basis, "rbf")){
    if (is.null(diff_basis$mus)){
      Mus <- matrix(NA_real_, nrow = N, ncol = diff_basis$M)
      colnames(Mus) <- paste("mu", seq(1, diff_basis$M), sep = "")
    }
  }

  diff_meth <- matrix(NA_real_, nrow = N, ncol = x_evals)
  xs_mat    <- matrix(NA_real_, nrow = N, ncol = x_evals)

  for (i in 1:N){
    # TODO: Choose more efficiently the x points or better?
    xs_mat[i, ] <- seq(from = min(contr_extr[i,1], treat_extr[i,1]),
                       to   = max(contr_extr[i,2], treat_extr[i,2]),
                       length.out = x_evals)

    if (is(contr_basis, "rbf")){
      if (is.null(contr_basis$mus)){
        contr_basis$mus <- contr_Mus[i, ]
        treat_basis$mus <- treat_Mus[i, ]
      }
    }

    # Evaluate the probit function for control samples
    y_contr <- eval_probit_function(contr_basis, xs_mat[i, ], contr_W[i, ])

    # Evaluate the probit function for treatment samples
    y_treat <- eval_probit_function(treat_basis, xs_mat[i, ], treat_W[i, ])

    # Obtain the difference bewteen the two methylation profiles
    diff_meth[i, ] <- y_contr - y_treat

    # Create the design matrix
    des_mat <- design_matrix(x = diff_basis, obs = xs_mat[i, ])

    # Create a data frame containing the regressors x and the target value y
    dataset <- data.frame(des_mat$H[, 2:NCOL(des_mat$H)], y = diff_meth[i, ])

    # Fit the model to the data using lm
    model <- lm(formula = y ~ ., data = dataset)

    if (is.null(fit_feature)){
      w <- coef(model)
    }else{
      w <- c(coef(model), summary(model)$sigma)
    }

    # Store the fitted parameters to the W_opt matrix
    W_opt[i, ] <- w

    if (is(diff_basis, "rbf")){
      if (is.null(diff_basis$mus)){
        Mus[i, ]   <- des_mat$basis$mus
      }
    }
  }
  return(list(W_opt      = W_opt,
              Mus        = Mus,
              xs_mat     = xs_mat,
              diff_meth  = diff_meth,
              diff_basis = diff_basis))
}
