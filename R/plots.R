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


#' Plot method for differential methylation
#'
#' Plot data using the differntial methylation data
#'
#' @param x Object of class 'diff_mpgex'
#' @param t Region to plot
#'
#' @export
plot.diff_mpgex <- function(x, t, ...){
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
