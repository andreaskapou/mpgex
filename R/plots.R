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
  # Draw a line evaluated from the polynomial function
  lines(x = xs, y=eval_probit_function(x, xs, w=w), col=2, lwd=2)
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
  # Draw a line evaluated from the rbf function
  lines(x=xs, y=eval_probit_function(x   = x,
                                     obs = xs,
                                     w   = w), col=2, lwd=2)
}

#' Plot method for differential methylation
#'
#' Plot data using the differntial methylation data
#'
#' @inheritParams plot.polynomial
#'
#' @export
plot.diff_meth <- function(x, obs, w, ...){
  assertthat::assert_that(is.matrix(obs))
  xmin <- min(obs[,1]) - 0.05
  xmax <- max(obs[,1]) + 0.05
  xs <- seq(xmin, xmax, len = 500)
  plot(obs[,1], obs[,3] / obs[,2],
       col = "darkgreen", pch = 24, xlim = c(xmin, xmax),
       xlab = "region x", ylab = "diff methylation level", ...)
  # Draw a line evaluated from the rbf function
  lines(x=xs, y=eval_function(x, xs, w), col=2, lwd=2)
}
