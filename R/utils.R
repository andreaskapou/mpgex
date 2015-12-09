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
