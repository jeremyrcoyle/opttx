#learner functions that might help to pick up effect modification, and therefore find a more optimal tx
#' @export
SL.glmem <- function(Y, X, newX, family, obsWeights, ...) {
  fit.glm <- glm(Y ~ A * ., data = X, family = family, weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
SL.glmnetem <- function(Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10, nlambda = 100, useMin = TRUE, ...) {
  require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + A * ., X)
    newX <- model.matrix(~-1 + A * ., newX)
  }
  fitCV <- cv.glmnet(x = X, y = Y, weights = obsWeights, lambda = NULL, type.measure = "deviance", nfolds = nfolds, 
                     family = family$family, alpha = alpha, nlambda = nlambda)
  pred <- predict(fitCV$glmnet.fit, newx = newX, s = ifelse(useMin, fitCV$lambda.min, fitCV$lambda.1se), type = "response")
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnetem"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' @export
predict.SL.glmnetem <- function(object, newdata, ...) {
  if (!is.matrix(newdata)) {
    newdata <- model.matrix(~-1 + A * ., newdata)
  }
  pred <- predict(object$object$glmnet.fit, newx = newdata, s = ifelse(object$useMin, object$object$lambda.min, object$object$lambda.1se), 
                  type = "response")
  return(pred)
}

