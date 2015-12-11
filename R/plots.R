#' S3 plot method for polynomial objects
#'
#' Plot data using the polynomial function as the response to observations x
#'
#' @param x Object type, which actually is the basis object.
#' @param dataset An L x 3 matrix containing the observations data.
#' @param w Vector of length M, containing the fitted parameters / coefficients
#'  for the M basis functions.
#' @param ... Additional parameters
#'
#' @export
plot.polynomial <- function(x, dataset, w, ...){
  assertthat::assert_that(is.matrix(dataset))
  xmin <- min(dataset[,1]) - 0.05
  xmax <- max(dataset[,1]) + 0.05
  xs <- seq(xmin, xmax, len = 500)
  plot(dataset[,1], dataset[,3] / dataset[,2],
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
#' @param mus The centers of the rbf function.
#'
#' @export
plot.rbf <- function(x, dataset, w, mus, ...){
  assertthat::assert_that(is.matrix(dataset))
  xmin <- min(dataset[,1]) - 0.05
  xmax <- max(dataset[,1]) + 0.05
  xs <- seq(xmin, xmax, len = 500)
  plot(dataset[,1], dataset[,3] / dataset[,2],
       col = "darkgreen", pch = 24, xlim = c(xmin, xmax), ylim = c(0,1),
       xlab = "region x", ylab = "methylation level")
  # Draw a line evaluated from the rbf function
  lines(x=xs, y=eval_probit_function(x = x,
                                     obs = xs,
                                     w = w,
                                     mus = mus), col=2, lwd=2)
}
