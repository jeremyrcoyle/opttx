
#' @export        
dV_from_preds <- function(preds) {
    max.col(preds)
}

# fit Q on A and subset of covariates V
learn_rule <- function(data, folds, nodes, split_preds, full_preds, predictions, 
    parallel = F, SL.library, verbose, ...) {
    
    if (length(nodes$Vnodes) == 1) {
        # glmnet expects at least two covariates, so don't use glmnet if we only have
        # one!
        SL.library$QaV <- setdiff(SL.library$QaV, "SL.glmnet")
        SL.library$class <- setdiff(SL.library$class, "mnSL.glmnet")
    } else if (length(nodes$Vnodes) == 0) {
        
        # if we have no covariates, there's no point in using learners that depend on
        # covariates
        SL.library$QaV <- c()
        SL.library$class <- "SL.mean"
    }
    
    QaV_fit <- NULL
    if (length(SL.library$QaV) > 0) {
        message_verbose("Fitting QaV", 2, verbose)
        
        QaV_fit <- origami_SuperLearner(folds = folds, data[, nodes$Ynode], data[, 
            nodes$Vnodes, drop = F], split_preds = split_preds, full_preds = full_preds, 
            SL.library = SL.library$QaV, family = gaussian(), cvfun = QaV_cv_SL, 
            .parallel = parallel, method = method.mvSL(method.NNLS()), ...)
        
        
    }
    
    class_fit <- NULL
    if (length(SL.library$class) > 0) {
        message_verbose("Fitting classfication", 2, verbose)
        class_fit <- multinomial_SuperLearner(folds = folds, data[, nodes$Ynode], 
            data[, nodes$Vnodes, drop = F], split_preds = split_preds, full_preds = full_preds, 
            SL.library = SL.library$class, cvfun = class_cv_SL, .parallel = parallel, 
            ...)
        
    }
    
    message_verbose("Fitting Combination (Refitting Weights)", 2, verbose)
    joint_fit <- joint_sl(QaV_fit, class_fit, predictions, data, nodes, risk_generator = create_tmle_risk)
    
    rule_object <- list(QaV_fit = QaV_fit, class_fit = class_fit, joint_fit = joint_fit, 
        nodes = nodes)
    class(rule_object) <- "opttx_rule"
    return(rule_object)
}


#' @export
predict.opttx_rule <- function(rule_object, newdata = "cv-original", pred_fit = c("joint", 
    "QaV", "class"), return_assignment = TRUE) {
    
    pred_fit <- match.arg(pred_fit)
    
    if (pred_fit == "QaV") {
        preds <- predict(rule_object$QaV_fit, newdata)$pred
    } else if (pred_fit == "class") {
        preds <- predict(rule_object$class_fit, newdata)$pred
    } else {
        preds <- predict(rule_object$joint_fit, rule_object$QaV_fit, rule_object$class_fit, 
            newdata)
    }
    
    if (return_assignment) {
        dV <- dV_from_preds(preds)
        return(dV)
    } else {
        return(preds)
    }
    
} 
