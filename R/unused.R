# currently unused functions, need to reincorporate into library

# function to construct blip dataset from fold-specific fits and then conduct superlearner to learn the blip function
# while respecting the CV V is vector of covariate names testing data fold=folds[[1]] Y=data$Y
# X=data[,c('A',W,missind)] family=binomial() SL.library=c('SL.glm','SL.glmnet','SL.rpart','SL.mean') id=1:nrow(X)
# obsWeights=1
blip_cv_SL <- function(fold, Y, X, SL.library, family, obsWeights, id, V, Q_fit, g_fit, use_full = F, ...) {
    v <- fold_index()
    train_idx <- training()
    valid_idx <- validation()
    if ((all.equal(train_idx, valid_idx) == T) || use_full) {
        # we're in final resubstitution call, so let's use full Q and g fits
        splitQ_fit <- Q_fit$fullFit
        splitg_fit <- g_fit$fullFit
    } else {
        # split specific Super Learners
        splitQ_fit <- Q_fit$foldFits[[v]]
        splitg_fit <- g_fit$foldFits[[v]]
    }
    # split specific estimates
    QAW <- predict(splitQ_fit, newdata = X)$pred
    new.data <- X
    new.data$A <- 0
    Q0W <- predict(splitQ_fit, newdata = new.data)$pred
    new.data$A <- 1
    Q1W <- predict(splitQ_fit, newdata = new.data)$pred
    
    # we don't want A in X from here on out glmnet gets picky if there's extra columns in the data
    A <- X$A
    X$A <- NULL
    pA1 <- predict(splitg_fit, X)$pred
    
    # split specific blip, class, and weights
    D1 <- (A/pA1 - (1 - A)/(1 - pA1)) * (Y - QAW) + Q1W - Q0W
    Z <- as.numeric(D1 < 0)
    K <- as.vector(abs(D1))  #D1 is a matrix somehow
    
    # should probably be more careful with obsWeights here fit split-specific blip based on split-specific Q and g
    cv_SL(fold, Z, X[, V, drop = F], SL.library, family, obsWeights * K, id, ...)
    # cv_SL(fold, D1, X[,V,drop=F], SL.library, family, obsWeights, id, ...)
}

# fully fits Q and g SLs in each fold of rule SL
nested_blip_cv_SL <- function(fold, Y, X, SL.library, family, obsWeights, id, V, fitQ, fitg, ...) {
    v <- fold_index()
    train_idx <- training()
    valid_idx <- validation()
    trainX <- training(X)
    trainY <- training(Y)
    AYstrata <- sprintf("%s %s", trainX$A, trainY)
    nestfolds <- make_folds(strata_ids = AYstrata, V = 10)
    
    
    splitQ_fit <- fitQ(nestfolds, trainY, trainX)
    splitg_fit <- fitg(nestfolds, trainX)
    # split specific estimates
    QAW <- predict(splitQ_fit, newdata = X)$pred
    new.data <- X
    new.data$A <- 0
    Q0W <- predict(splitQ_fit, newdata = new.data)$pred
    new.data$A <- 1
    Q1W <- predict(splitQ_fit, newdata = new.data)$pred
    
    # we don't want A in X from here on out glmnet gets picky if there's extra columns in the data
    A <- X$A
    X$A <- NULL
    pA1 <- predict(splitg_fit, X)$pred
    
    # split specific blip, class, and weights
    D1 <- (A/pA1 - (1 - A)/(1 - pA1)) * (Y - QAW) + Q1W - Q0W
    Z <- as.numeric(D1 < 0)
    K <- as.vector(abs(D1))  #D1 is a matrix somehow
    
    # should probably be more careful with obsWeights here fit split-specific blip based on split-specific Q and g
    cv_SL(fold, Z, X[, V, drop = F], SL.library, family, obsWeights * K, id, ...)
    # cv_SL(fold, D1, X[,V,drop=F], SL.library, family, obsWeights, id, ...)
} 
