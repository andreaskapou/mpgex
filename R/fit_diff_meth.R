#' Fit differential methylation profiles using BLM
#'
#' Function for fitting the differential methylation profiles and getting back
#' the fitted coefficients using the Basis Linear Model.
#'
#' @param control The optimized values for the control samples.
#' @param treatment The optimized values for the treatment samples.
#' @param diff_basis The basis object for differential methylation.
#' @param lambda The regularization parameter when fitting BLM.
#'
#' @return The fitted coefficients
#'
#' @export
fit_diff_meth <- function(control, treatment, diff_basis, lambda){
  N <- NROW(control$W_opt)
  W_opt <- matrix(NA_real_, nrow = N, ncol = diff_basis$M +1)
  colnames(W_opt) <- paste("w", seq(1, diff_basis$M + 1), sep = "")

  for (i in 1:N){
    # TODO: Choose more efficiently the x points
    xs <- seq(min(control$x_extrema[i,1], treatment$x_extrema[i,1]),
              max(control$x_extrema[i,2], treatment$x_extrema[i,2]),
              length.out = 50)

    if (is(control$basis, "rbf")){
      control$basis$mus   <- control$Mus[i, ]
      treatment$basis$mus <- treatment$Mus[i, ]
    }
    # Evaluate the probit function
    y_contr <- eval_probit_function(control$basis, xs, control$W_opt[i, ])
    y_treat <- eval_probit_function(treatment$basis, xs, treatment$W_opt[i, ])
    # Obtain the difference bewteen the two methylation profiles
    diff_meth <- y_contr - y_treat
    # Fit the model to the data using the Basis Linear Model
    model_fit <- blm(x = xs,
                     y = diff_meth,
                     basis  = diff_basis,
                     lambda = lambda,
                     return.all = FALSE)
    W_opt[i, ] <- coef(model_fit)
  }
  return(list(W_opt = W_opt))
}
