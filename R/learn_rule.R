
# fit on subset of covariates V
learn_rule <- function(data, folds, nodes, split_preds, SL.library = c("SL.glm", "SL.glmnet", "SL.step.forward", "SL.gam", 
    "SL.rpart", "SL.rpartPrune", "SL.mean"), blip_type, maximize = T, parallel = F) {
    
    if (length(nodes$Vnodes) == 1) {
        # glmnet expects at least two covariates, so don't use glmnet if we only have one!
        SL.library <- setdiff(SL.library, "SL.glmnet")
    } else if (length(nodes$Vnodes) == 0) {
        # if we have no covariates, there's no point in using learners that depend on covariates
        SL.library <- "SL.mean"
    }
    
    if (blip_type == "QbV.mse") {
        method <- method.NNLS()
        family <- gaussian()
    } else {
        method <- method.surlog()
        family <- binomial()
    }
    
    blip_fit <- origami_SuperLearner(folds = folds, data[, nodes$Ynode], data[, nodes$Vnodes, drop = F], split_preds = split_preds, 
        SL.library = SL.library, family = family, cvfun = split_cv_SL, method = method, cts.num = 10, , .parallel = parallel, 
        blip_type = blip_type, maximize = maximize)
    
    return(blip_fit)
} 
