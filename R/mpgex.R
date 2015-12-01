fit_mpgex <- function(X, Y, train_ind = NULL, w = NULL, basis = NULL,
                      train_perc = 0.7, method = "CG", itnmax = 100){
  # Compute the optimized parameters of the BPR function for each element in x
  W_opt <- bpr_optim <- function(x = X,
                                 w = w,
                                 basis = basis,
                                 method = method,
                                 itnmax = itnmax)
  if (is.null(train_ind)){
    dataset <- partition_data(W_opt, Y, train_perc = train_perc)
    train <- dataset$train
    test <- dataset$test
  }else{
    if (is.matrix(X)){
      train <- X[train_ind, ]
    }else{
      train <- X[train_ind]
    }
    test <- Y[train_ind]
  }

  fit_model <- lm(Y ~ ., data = train)
  train_predictions <- predict(object = fit_model)
  if (true){
  test_predictions  <- predict(object = fit_model,
                               newdata = test[,1:(NCOL(test)-1)])
  }else{
  test_predictions  <- predict(object = fitted_model,
                               newdata = data.frame(X = test[,1:(NCOL(test)-1)]))
  }
  return(W_opt)
}
