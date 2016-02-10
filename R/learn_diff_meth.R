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
#' @param x_eval Integer denoting the number of evaluation points.
#' @param is_parallel Logical, indicating if code should be run in parallel.
#' @param no_cores Number of cores to be used, default is max_no_cores - 1.
#'
#' @return The differential methylation profiles and the corresponding fitted
#'  coefficients from the BLM model.
#'
#' @importFrom methods is
#' @export
learn_diff_meth <- function(control, treatment, diff_basis, fit_feature = NULL,
                            lambda = 0, x_eval = 50, is_parallel = is_parallel,
                            no_cores = no_cores){
  # Create matrices from lists for efficiency
  contr_W     <- control$W_opt
  treat_W     <- treatment$W_opt
  contr_extr  <- control$x_extrema
  treat_extr  <- treatment$x_extrema

  contr_basis <- control$basis
  treat_basis <- treatment$basis

  contr_Mus <- NULL
  treat_Mus <- NULL
  # If we have 'rbf' basis, extract the Mus parameters
  if (methods::is(contr_basis, "rbf")){
    if (is.null(contr_basis$mus)){
      contr_Mus <- control$Mus
      treat_Mus <- treatment$Mus
    }
  }
  # Extract total number of observations
  N <- NROW(contr_W)
  # Initialize i so the CMD check for R does not complain
  i <- 0

  # If parallel mode is ON
  if (is_parallel){
    # If number of cores is not given
    if (is.null(no_cores)){
      no_cores <- parallel::detectCores() - 1
    }
    if (is.na(no_cores)){
      no_cores <- 2
    }
    # Create cluster object
    cl <- parallel::makeCluster(no_cores)
    doParallel::registerDoParallel(cl)

    # Parallel fitting for each element  i.e. for each region i.
    res <- foreach::"%dopar%"(obj = foreach::foreach(i = 1:N),
                              ex  = {
                        out <- .fit_diff_meth(contr_W     = contr_W[i, ],
                                              treat_W     = treat_W[i, ],
                                              contr_basis = contr_basis,
                                              treat_basis = treat_basis,
                                              diff_basis  = diff_basis,
                                              fit_feature = fit_feature,
                                              contr_extr  = contr_extr[i, ],
                                              treat_extr  = treat_extr[i, ],
                                              contr_Mus   = contr_Mus[i, ],
                                              treat_Mus   = treat_Mus[i, ],
                                              x_eval      = x_eval)
                              })
    # Stop parallel execution
    parallel::stopCluster(cl)
  }else{
    # Sequential fitting for each element i.e. for each region i.
    res <- foreach::"%do%"(obj = foreach::foreach(i = 1:N),
                            ex  = {
                      out <- .fit_diff_meth(contr_W     = contr_W[i, ],
                                            treat_W     = treat_W[i, ],
                                            contr_basis = contr_basis,
                                            treat_basis = treat_basis,
                                            diff_basis  = diff_basis,
                                            fit_feature = fit_feature,
                                            contr_extr  = contr_extr[i, ],
                                            treat_extr  = treat_extr[i, ],
                                            contr_Mus   = contr_Mus[i, ],
                                            treat_Mus   = treat_Mus[i, ],
                                            x_eval      = x_eval)
                              })
  }

  # Matrix for storing optimized coefficients
  W_opt <- sapply(res, function(x) x$W_opt)
  if (is.matrix(W_opt)){
    W_opt <- t(W_opt)
  }else{
    W_opt <- as.matrix(W_opt)
  }
  colnames(W_opt) <- paste("w", seq(1, NCOL(W_opt)), sep = "")

  # Matrix for storing the centers of RBFs if object class is 'rbf'
  Mus <- NULL
  if (methods::is(diff_basis, "rbf")){
    if (is.null(diff_basis$mus)){
      Mus <- sapply(res, function(x) x$Mus)
      if (is.matrix(Mus)){
        Mus <- t(Mus)
      }else{
        Mus <- as.matrix(Mus)
      }
      colnames(Mus) <- paste("mu", seq(1, NCOL(Mus)), sep = "")
    }
  }

  # Matrix for storing differential methylation values
  diff_meth <- sapply(res, function(x) x$diff_meth)
  if (is.matrix(diff_meth)){
    diff_meth <- t(diff_meth)
  }else{
    diff_meth <- as.matrix(diff_meth)
  }
  colnames(diff_meth) <- paste("y", seq(1, NCOL(diff_meth)), sep = "")


  # Matrix for storing x's values
  xs_mat <- sapply(res, function(x) x$xs_mat)
  if (is.matrix(xs_mat)){
    xs_mat <- t(xs_mat)
  }else{
    xs_mat <- as.matrix(xs_mat)
  }
  colnames(xs_mat) <- paste("x", seq(1, NCOL(xs_mat)), sep = "")

  return(list(W_opt      = W_opt,
              Mus        = Mus,
              xs_mat     = xs_mat,
              diff_meth  = diff_meth,
              diff_basis = diff_basis))
}


# Private function for fitting differential methylation profiles
# for a given region
.fit_diff_meth <- function(contr_W, treat_W, contr_basis, treat_basis,
                           diff_basis, fit_feature, contr_extr, treat_extr,
                           contr_Mus, treat_Mus, x_eval){
  # TODO: Choose more efficiently the x points or better?
  xs_mat <- seq(from = min(contr_extr[1], treat_extr[1]),
                to   = max(contr_extr[2], treat_extr[2]),
                length.out = x_eval)

  if (methods::is(contr_basis, "rbf")){
    if (is.null(contr_basis$mus)){
      contr_basis$mus <- contr_Mus
      treat_basis$mus <- treat_Mus
    }
  }

  # Evaluate the probit function for control samples
  y_contr <- eval_probit_function(contr_basis, xs_mat, contr_W)

  # Evaluate the probit function for treatment samples
  y_treat <- eval_probit_function(treat_basis, xs_mat, treat_W)

  # Obtain the difference bewteen the two methylation profiles
  diff_meth <- y_contr - y_treat

  # Create the design matrix
  des_mat <- design_matrix(x = diff_basis, obs = xs_mat)

  # Create a data frame containing the regressors x and the target value y
  dataset <- data.frame(des_mat$H[, 2:NCOL(des_mat$H)], y = diff_meth)

  # Fit the model to the data using lm
  model <- stats::lm(formula = y ~ ., data = dataset)

  if (is.null(fit_feature)){
    w <- stats::coef(model)
  }else{
    w <- c(stats::coef(model), summary(model)$sigma)
  }

  Mus <- NULL
  if (methods::is(diff_basis, "rbf")){
    if (is.null(diff_basis$mus)){
      Mus <- des_mat$basis$mus
    }
  }
  return(list(W_opt = w,
              Mus = Mus,
              diff_meth = diff_meth,
              xs_mat = xs_mat))
}
