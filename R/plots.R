#' S3 plot method for polynomial objects
#'
#' Plot data using the polynomial function as the response to observations x
#'
#' @param x The basis object.
#' @param obs An L x 3 matrix containing the observations data.
#' @param w Vector of length M, containing the fitted parameters / coefficients
#'  for the M basis functions.
#' @param ... Additional parameters
#'
#' @export
plot.polynomial <- function(x, obs, w, ...){
  assertthat::assert_that(is.matrix(obs))
  xmin <- min(obs[,1]) - 0.05
  xmax <- max(obs[,1]) + 0.05
  xs <- seq(xmin, xmax, len = 500)
  plot(obs[,1], obs[,3] / obs[,2],
       col = "darkgreen", pch = 24, xlim = c(xmin, xmax), ylim = c(0,1),
       xlab = "region x", ylab = "methylation level")
  # Draw a curve evaluated from the basis function
  lines(x = xs, y = eval_probit_function(x, xs, w), col = 2, lwd = 2)
}


#' S3 plot method for rbf objects
#'
#' Plot data using the rbf function as the response to observations x
#'
#' @inheritParams plot.polynomial
#'
#' @export
plot.rbf <- function(x, obs, w, ...){
  assertthat::assert_that(is.matrix(obs))
  xmin <- min(obs[,1]) - 0.05
  xmax <- max(obs[,1]) + 0.05
  xs <- seq(xmin, xmax, len = 500)
  plot(obs[,1], obs[,3] / obs[,2],
       col = "darkgreen", pch = 24, xlim = c(xmin, xmax), ylim = c(0,1),
       xlab = "region x", ylab = "methylation level")
  # Draw a curve evaluated from the basis function
  lines(x = xs, y = eval_probit_function(x, xs, w), col = 2, lwd = 2)
}


#' Plot fit of methylation profiles across a region
#'
#' @param region Promoter region number
#' @param X Methylation data observations
#' @param fit_prof Fitted profiles
#' @param fit_mean Fitted mean function
#' @param title Title of the plot
#' @param ... Additional parameters
#'
#' @export
plot_fitted_profiles <- function(region, X, fit_prof, fit_mean = NULL,
                                 title = "Gene promoter", ...){

  par(cex=1.05, mai=c(1.37,1.37,.7,.3) )
  x <- X[[region]][,1]
  y <- X[[region]][,3]/X[[region]][,2]
  xs <- seq(from = -1, to = 1, by = 0.01)
  plot(x, y, col = "blue2", pch = 21, ylim = c(0,1), xlim = c(-1,1), lwd = 0.8,
       xlab = NA, ylab = NA, cex.axis = 1.1, xaxt = "n")
  mtext(side = 1, "genomic region", line = 3, cex = 1.2)
  mtext(side = 2, "methylation level", line = 3, cex = 1.2)
  axis(side = 1, at = c(-1, 0, 1), labels=c("-7kb", "TSS", "+7kb"))
  title(main=title, line = 1, cex.main=1.4)
  if(!is.null(fit_mean)){
    lines(x = xs,
          y = eval_probit_function(fit_mean$basis, xs, fit_mean$W_opt[region, ]),
          col = 'coral', lwd = 2, lty = 2)
  }
  lines(x = xs,
        y = eval_probit_function(fit_prof$basis, xs, fit_prof$W_opt[region, ]),
        col = 'red2', lwd = 2)
}

#' Plot method for differential methylation
#'
#' Plot data using the differntial methylation data
#'
#' @param x Object of class 'diff_mpgex'
#' @param t Region to plot
#' @param ... Additional parameters
#'
#' @export
plot.mpgex_diff_regr <- function(x, t, ...){
  basis <- x$diff_basis
  obs <- x$diff_meth[t, ]
  w <- x$W_opt[t, ]

  xmin <- min(obs) - 0.05
  xmax <- max(obs) + 0.05
  xs <- seq(xmin, xmax, len = 500)
  plot(obs, col = "darkgreen", pch = 24, xlim = c(xmin, xmax),
       xlab = "region x", ylab = "diff methylation level", ...)
  # Draw a curve evaluated from the basis function
  lines(x = xs, y = eval_function(basis, xs, w), col = 2, lwd = 2)
}


#' Plot method for differential methylation
#'
#' Plot data using the differntial methylation data
#'
#' @param x Object of class 'diff_mpgex'
#' @param t Region to plot
#' @param ... Additional parameters
#'
#'
#' @export
plot_diff_mpgex <- function(x, t, ...){
  basis <- x$diff_basis
  xs <- x$xs_mat[t, ]
  obs <- x$diff_meth[t, ]
  w <- x$W_opt[t, ]

  xmin <- min(xs) - 0.05
  xmax <- max(xs) + 0.05
  plot(xs, obs, col = "darkgreen", pch = 24, xlim = c(xmin, xmax),
       xlab = "region x", ylab = "diff methylation level", ...)
  # Draw a curve evaluated from the basis function
  xs <- seq(xmin, xmax, len = 500)
  lines(x = xs, y = eval_function(basis, xs, w), col = 2, lwd = 2)
}
