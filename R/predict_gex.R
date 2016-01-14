#' Predict gene expression data from methylation profiles
#'
#' \code{predict_gex} predicts gene expression data from methylation profiles,
#'  using linear regression. The same method can be used to predict
#'  differential expression from differential methylation profiles.
#'
#' @param formula An object of class \code{\link[stats]{formula}} needed when
#'  calling the \code{\link[stats]{lm}} function for performing linear
#'  regression. If NULL, the simple linear regression model is used.
#' @param model_name A charcter denoting the regression model
#' @param train The training data.
#' @param test The test data.
#' @param is_summary Logical, print the summary statistics.
#'
#' @return A list containing the following elements:
#' \itemize{
#'  \item \code{formula} the formula that was used.
#'  \item \code{gex_model} the fitted model.
#'  \item \code{train_pred} the training predicted values.
#'  \item \code{test_pred} the test predicted values.
#'  \item \code{train_errors} the training error metrics.
#'  \item \code{test_errors} the test error metrics.
#' }
#'
#' @seealso \code{\link{calculate_errors}}
#'
#' @export
predict_gex <- function(formula = NULL, model_name = "lm", train, test, is_summary = TRUE){
  if (is.null(formula)){
    formula <- Y ~ .
  }
  if (model_name == "mars"){
    model <- earth::earth(formula = formula, data = train)
  }else if (model_name == "randomForest"){
    model <- randomForest::randomForest(formula = formula, data = train)
  }else if (model_name == "rlm"){
    model <- MASS::rlm(formula = formula, data = train)
  }else if (model_name == "svm"){
    model <- e1071::svm(formula = formula, data = train, kernel='radial', cross=10)
  }else{
    model <- lm(formula = formula, data = train)
  }
  #train_pred <- fitted(object = model)
  train_pred <- predict(object = model,
                        type = "response")

  regressors <- 1:(NCOL(test) - 1)
  if (! NCOL(test) == 2){
    test_pred <- predict(object  = model,
                         newdata = test[ ,regressors],
                         type = "response")
  }else{
    test_pred <- predict(object  = model,
                         newdata = data.frame(X = test[ ,regressors]),
                         type = "response")
  }

  # Calculate model errors
  message("-- Train Errors --")
  train_errors <- calculate_errors(train$Y, train_pred, summary = is_summary)
  message("-- Test Errors --")
  test_errors  <- calculate_errors(test$Y,  test_pred,  summary = is_summary)

  out <- list(formula      = formula,
              gex_model    = model,
              train_pred   = train_pred,
              test_pred    = test_pred,
              train_errors = train_errors,
              test_errors  = test_errors)
  return(out)
}
