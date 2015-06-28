
# split specific predictions for Q and g and derived quantities
cv_split_preds <- function(fold, data, nodes, fits, use_full = F, maximize = T) {
  v <- fold_index()
  train_idx <- training()
  valid_idx <- validation()
  
  if ((all.equal(train_idx, valid_idx) == T) || use_full) {
    # we're in final resubstitution call, so let's use full Q and g fits
    splitQ_fit <- fits$Q_fit$fullFit
    splitg_fit <- fits$g_fit$fullFit
  } else {
    # split specific Super Learners
    splitQ_fit <- fits$Q_fit$foldFits[[v]]
    splitg_fit <- fits$g_fit$foldFits[[v]]
  }
  
  # split specific estimates
  new_data <- data[,c(nodes$Anode,nodes$Wnodes)]
  QAW <- predict(splitQ_fit, newdata = new_data)$pred
  new_data[,nodes$Anode] <- 0
  Q0W <- predict(splitQ_fit, newdata = new_data)$pred
  new_data[,nodes$Anode] <- 1
  Q1W <- predict(splitQ_fit, newdata = new_data)$pred
  
  pA1 <- predict(splitg_fit, new_data[,nodes$Wnodes])$pred
  
  # split specific blip, class, and weights
  A=data[,nodes$Anode]
  Y=data[,nodes$Ynode]
  D1 <- (A/pA1 - (1 - A)/(1 - pA1)) * (Y - QAW) + Q1W - Q0W
  if(maximize){
    Z <- as.numeric(D1 > 0)
  } else {
    Z <- as.numeric(D1 < 0)
  }
  K <- as.vector(abs(D1))  #D1 is a matrix somehow
  
  list(QAW=QAW, Q1W=Q1W, Q0W=Q0W, pA1=pA1, D1 = D1, Z = Z, K = K)
}

#ignores Y
#and instead uses Z generated from a split-specific training set
# generated using cv_split_preds
# X are the nodes to base the rule on
split_cv_SL <- function(fold, Y, X, SL.library, family, obsWeights, id, use_full = F, split_preds, ...) {
  v <- fold_index()
  train_idx <- training()
  valid_idx <- validation()
  
  # should probably be more careful with obsWeights here 
  # fit split-specific blip based on split-specific Q and g
  cv_SL(fold, split_preds$Z[[v]], X, SL.library, family, obsWeights * split_preds$K[[v]], id, ...)
}


# alex's log loss for weighted classification based opt tx approach
method.surlog <- function() {
  out <- list(require = NULL, computeCoef = function(Z, Y, libraryNames, verbose, obsWeights, ...) {
    cvRisk <- apply(Z, 2, function(x) {
      -mean(obsWeights * ifelse(Y, log(x), log(1 - x)))
    })
    names(cvRisk) <- libraryNames
    
    preds <- Z - 1/2
    wgts <- obsWeights
    risk.fun <- function(b) {
      mean(wgts * (-plogis((2 * Y - 1) * (c(preds %*% cbind(b))), log.p = TRUE)))
    }
    num.alg <- ncol(preds)
    if (num.alg == 1) {
      coef <- 1
    } else {
      init <- rep(1/num.alg, num.alg)
      alpha.out <- optim(init, risk.fun, method = "L-BFGS-B", lower = rep(0, num.alg), upper = rep(1, num.alg))$par
      coef <- (alpha.out/sum(alpha.out))
    }
    out <- list(cvRisk = cvRisk, coef = coef)
    return(out)
  }, computePred = function(predY, coef, ...) {
    out <- predY %*% coef
    return(out)
  })
  invisible(out)
}

