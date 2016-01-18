#' Compute the min-max scaling
#'
#' \code{minmax_scaling} normalizes a given vector using the the min-max
#' scaling method. More formally:
#' \deqn{scaled = \frac{data -x_{min}}{x_{max} - x_{min}} \times (f_{max} - f_{min}) + f_{min}}
#'
#' @param data Vector with numeric data to be scaled.
#' @param xmin Optional minimum value, otherwise \code{min(data)} will be used.
#' @param xmax Optional maximum value, otherwise \code{max(data)} will be used.
#' @param fmin Optional minimum range value, default is 0.
#' @param fmax Optional maximum range value, default is 1.
#'
#' @return The scaled data in the given range, default is between (0, 1). If
#'  xmin = xmax the input vector \code{data} is returned.
#'
#' @examples
#' data <- c(-20, 0, 15, 20)
#' scaled <- minmax_scaling(data)
#'
#' @export
minmax_scaling <- function(data, xmin = NULL, xmax = NULL, fmin = 0, fmax = 1){
  if (is.null(xmin)){
    xmin <- min(data)
  }
  if (is.null(xmax)){
    xmax <- max(data)
  }
  if ( (xmin - xmax) == 0){
    return(data)
  }
  minmax <- (data - xmin) / (xmax - xmin)
  minmax_scaled <- minmax * (fmax - fmin) + fmin

  return(minmax_scaled)
}

#' Partition data in train and test set
#'
#' \code{partition_data} partition data randomly in train and test sets.
#'
#' @param X Input / Independent variables
#' @param Y Dependent variables
#' @param train_ind Index of traininig data, if NULL a random one is generated.
#' @param train_perc Percentage of training data when partitioning.
#'
#' @return A list containing the train, test data and index of training data.
#'
#' @examples
#' X <- matrix(c(-20, 0, 15, 20, 20, 10, 4, 5), ncol=2)
#' Y <- c(2, 5, 0, 8)
#' data <- partition_data(X, Y)
#'
#' ##
#' X <- matrix(c(-20, 0, 15, 20, 20, 10, 4, 5), ncol=2)
#' Y <- c(2, 5, 0, 8)
#' train_ind <- c(1,4)
#' data <- partition_data(X, Y, train_ind)
#'
#' @export
partition_data <- function(X, Y, train_ind = NULL, train_perc = 0.7){
  if (is.null(train_ind)){
    pivot <- NROW(X) * train_perc
    train_ind <- sample(NROW(X), round(pivot))
  }

  if (is.matrix(X)){
    train <- data.frame(X = X[train_ind, ],
                        Y = as.matrix(Y[train_ind]))
    test  <- data.frame(X = X[-train_ind, ],
                        Y = as.matrix(Y[-train_ind]))
  }else{
    train <- data.frame(X = X[train_ind],
                        Y = as.matrix(Y[train_ind]))
    test  <- data.frame(X = X[-train_ind],
                        Y = as.matrix(Y[-train_ind]))
  }
  return(list(train = train, test = test, train_ind = train_ind))
}


#' Calculate error metrics
#'
#' \code{calculate_errors} calculates error metrics so as to assess the
#'  performance of a model, e.g. linear regression.
#' @param x Actual values.
#' @param y Predicted values.
#' @param summary Logical, indicating if the erorr metrics should be printed.
#'
#' @return A list containing the following components:
#' \itemize{
#'  \item \code{mae} mean absolute error.
#'  \item \code{mse}  mean squared error (the variance).
#'  \item \code{rmse} root mean squared error (std dev).
#'  \item \code{mape} mean absolute percentage error.
#'  \item \code{rstd} relative standard deviation.
#' }
#'
#' @examples
#' # Example data
#' actual <- c(4, 6, 9, 10, 4, 6, 4, 7, 8, 7)
#' predicted <- c(5, 6, 8, 10, 4, 8, 4, 9, 8, 9)
#'
#' calculate_errors(actual, predicted, summary = TRUE)
#'
#' @export
calculate_errors <- function(x, y, summary = FALSE){
  R <- list()
  if (! is.numeric(x) || ! is.numeric(y))
    stop("Arguments 'x' and 'y' must be numeric vectors.")
  if (length(x) != length(y))
    stop("Vectors 'x' and ' y' must have the same length.")
  error  <- y - x
  # mean absolute error
  R$mae  <- mean(abs(error))
  mae_f  <- formatC(R$mae, digits = 4, format = "f")
  # mean squared error (the variance?!)
  R$mse  <- mean(error ^ 2)
  mse_f  <- formatC(R$mse, digits = 4, format = "f")
  # root mean squared error (std. dev.)
  R$rmse <- sqrt(R$mse)
  rmse_f <- formatC(R$rmse, digits = 4, format = "f")
  # mean absolute percentage error
  R$mape <- mean(abs(error / x))
  mape_f <- formatC(R$mape, digits = 4, format = "f")
  # relative standard deviation
  R$rstd <- R$rmse / mean(x)
  rstd_f <- formatC(R$rstd, digits = 4, format = "f")
  # R-Squared
  R$rsq  <- 1 - (sum(error ^ 2) / sum((x - mean(x)) ^ 2))
  rsq_f  <- formatC(R$rsq, digits = 4, format = "f")
  # Pearson Correlation Coefficient
  R$pcc  <- cor(x, y)
  pcc_f  <- formatC(R$pcc, digits = 4, format = "f")
  if (summary) {
    cat("-- Error Terms ----\n")
    cat(" MAE:  ", mae_f, "\n")
    cat(" MSE:  ", mse_f, "\n")
    cat(" RMSE: ", rmse_f, "\n")
    cat(" MAPE: ", mape_f, "\n")
    cat(" rSTD: ", rstd_f, "\n")
    cat(" R-sq: ", rsq_f, "\n")
    cat(" PCC:  ", pcc_f, "\n")
    cat("-------------------\n\n")
  }
  if (summary) {
    invisible(R)
  }else{
    return(R)
  }
}
