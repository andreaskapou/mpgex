#' Fit differential methylation profiles using BLM
#'
#' Function for fitting the differential methylation profiles and returning
#' the fitted coefficients using the Basis Linear Model.
#'
#' @param control The optimized values for the control samples.
#' @param treatment The optimized values for the treatment samples.
#' @param diff_basis The basis object for differential methylation.
#' @param lambda The regularization parameter when fitting the BLM.
#' @param evals Integer denoting the number of evaluation points.
#'
#' @return The differential methylation profiles and the corresponding fitted
#'  coefficients from the BLM.
#'
#' @export
fit_diff_meth <- function(control, treatment, diff_basis, lambda, evals = 50){

  # Create matrices from lists for efficiency
  contr_W     <- control$W_opt
  treat_W     <- treatment$W_opt
  contr_extr  <- control$x_extrema
  treat_extr  <- treatment$x_extrema
  contr_basis <- control$basis
  treat_basis <- treatment$basis
  if (is(control$basis, "rbf")){
    contr_Mus <- control$Mus
    treat_Mus <- treatment$Mus
  }

  N <- NROW(control$W_opt)
  W_opt <- matrix(NA_real_, nrow = N, ncol = diff_basis$M + 1)
  colnames(W_opt) <- paste("w", seq(1, diff_basis$M + 1), sep = "")

  diff_meth <- matrix(NA_real_, nrow = N, ncol = evals)
  xs_mat    <- matrix(NA_real_, nrow = N, ncol = evals)

  for (i in 1:N){
    # TODO: Choose more efficiently the x points
    xs_mat[i, ] <- seq(min(contr_extr[i,1], treat_extr[i,1]),
                       max(contr_extr[i,2], treat_extr[i,2]),
                       length.out = evals)

    if (is(control$basis, "rbf")){
      contr_basis$mus <- contr_Mus[i, ]
      treat_basis$mus <- treat_Mus[i, ]
    }
    # Evaluate the probit function
    y_contr <- eval_probit_function(contr_basis, xs_mat[i, ], contr_W[i, ])
    y_treat <- eval_probit_function(treat_basis, xs_mat[i, ], treat_W[i, ])
    # Obtain the difference bewteen the two methylation profiles
    diff_meth[i, ] <- y_contr - y_treat
    # Fit the model to the data using the Basis Linear Model
    model_fit <- blm(x = xs_mat[i, ],
                     y = diff_meth[i, ],
                     basis  = diff_basis,
                     lambda = lambda,
                     return.all = FALSE)
    W_opt[i, ] <- coef(model_fit)
  }
  return(list(W_opt = W_opt,
              xs_mat = xs_mat,
              diff_meth = diff_meth,
              diff_basis = diff_basis))
}
