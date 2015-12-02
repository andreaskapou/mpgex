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
#' @param train_perc Percentage of training data when partitioning.
#'
#' @return A list containing the train and test data.
#'
#' @examples
#' X <- matrix(c(-20, 0, 15, 20), ncol=2)
#' Y <- c(2,5)
#' data <- partition_data(X, Y)
#'
#' @export
partition_data <- function(X, Y, train_perc = 0.7){
  pivot <- NROW(X) * train_perc
  shuffle <- sample(NROW(X))
  if (is.matrix(X)){
    X <- X[shuffle, ]
  }else{
    X <- X[shuffle]
  }
  Y <- Y[shuffle]

  if (is.matrix(X)){
    train <- data.frame(X = X[1:pivot, ],
                        Y = as.matrix(Y[1:pivot]))
    test  <- data.frame(X = X[(pivot + 1):NROW(X), ],
                        Y = as.matrix(Y[(pivot + 1):NROW(X)]))
  }else{
    train <- data.frame(X = X[1:pivot],
                        Y = as.matrix(Y[1:pivot]))
    test  <- data.frame(X = X[(pivot + 1):NROW(X)],
                        Y = as.matrix(Y[(pivot + 1):NROW(X)]))
  }

  return(list(train = train, test = test))
}



#'
#' Polynomial function of degree M.
#' X is the input vector and w are the coefficients of the
#' polynomial whose length is M + 1
#'
#'@export
eval_prob_polyn_func = function(x, w){
  f <- rep(0, length(x))
  M <- length(w)
  for (i in 1:M){
    f <- f + w[i] * polynomial_basis(x, i - 1)
  }
  Phi <- pnorm(f)
  return(Phi)
}


#'
#' Radia basis function of degree M.
#' X is the observation
#'
#'@export
eval_prob_rbf_func = function(x, w, basis, mus){
  f <- rep(w[1], NROW(x))
  M <- length(mus)
  x <- as.matrix(x)
  for (i in 1:M){
    f <- f + w[i + 1] * apply(x, 1, rbf_basis, mus = mus[i], gamma = basis$gamma)
  }
  Phi <- pnorm(f)
  return(Phi)
}

#'
#' Plot function for polynomial
#'
#' @export
plot.polynomial <- function(basis, dataset, w){
  xmin <- min(dataset[,1]) - 0.05
  xmax <- max(dataset[,1]) + 0.05
  xs <- seq(xmin, xmax, len = 500)
  plot(dataset[,1], dataset[,3] / dataset[,2],
       col = "darkgreen", pch = 24, xlim = c(xmin, xmax), ylim = c(0,1),
       xlab = "region x", ylab = "methylation level")
  lines(x = xs, y=eval_prob_polyn_func(xs, w=w), col=2, lwd=2)
}


#'
#' Plot function for rbf
#'
#' @export
plot.rbf <- function(basis, dataset, w, mus){
  xmin <- min(dataset[,1]) - 0.05
  xmax <- max(dataset[,1]) + 0.05
  xs <- seq(xmin, xmax, len = 500)
  plot(dataset[,1], dataset[,3] / dataset[,2],
       col = "darkgreen", pch = 24, xlim = c(xmin, xmax), ylim = c(0,1),
       xlab = "region x", ylab = "methylation level")
  lines(x=xs, y=eval_prob_rbf_func(xs, w = w, basis = basis, mus = mus), col=2, lwd=2)
}
