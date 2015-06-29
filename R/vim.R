risk_from_osl <- function(osl_fit) {
    function(pred, Y) {
        method <- osl_fit$fullFit$method
        coef <- method$computeCoef(pred, Y, libraryNames = "split_SL", obsWeights = osl_fit$valWeights, verbose = F)
        risk <- coef$cvRisk
    }
}


#' @export
#'
SL_vim <- function(osl_fit, data, Xnodes, ...) {
    full_preds <- cv_predict_original(osl_fit)
    risk_fun <- risk_from_osl(osl_fit)
    Y <- osl_fit$valY
    full_risk <- risk_fun(full_preds, Y)
    riskdf <- ldply(seq_along(Xnodes), .progress = "text", function(idx, ...) {
        new_Xnodes <- Xnodes[-1 * idx]
        new_fit <- origami_SuperLearner(Y = Y, X = data[, new_Xnodes], folds = osl_fit$folds, family = osl_fit$fullFit$family, 
            SL.library = osl_fit$SL.library, ...)
        new_preds <- cv_predict_original(new_fit)
        risk_from_Y <- risk_fun(new_preds, Y)
        risk_from_full <- risk_fun(new_preds, full_preds)
        p_disagree <- mean((full_preds > 0.5) != (new_preds > 0.5))
        data.frame(node = Xnodes[idx], risk_from_Y = risk_from_Y, p_disagree = p_disagree, risk_from_full = risk_from_full)
    })
    
    no_fit <- origami_SuperLearner(Y = Y, X = data[, Xnodes], folds = osl_fit$folds, family = osl_fit$fullFit$family, 
        SL.library = "SL.mean", ...)
    no_pred <- cv_predict_original(no_fit)
    no_info_Y <- risk_fun(no_pred, Y)
    no_p_disagree <- mean((full_preds > 0.5) != (no_pred > 0.5))
    no_info_full <- risk_fun(no_pred, full_preds)
    no_info <- data.frame(node = "empty", risk_from_Y = no_info_Y, p_disagree = no_p_disagree, risk_from_full = no_info_full)
    riskdf <- rbind(riskdf, no_info)
    riskdf$risk_percent <- (riskdf$risk_from_Y - full_risk)/riskdf$risk_from_Y
    riskdf$risk_full_fraction <- (riskdf$risk_from_full)/no_info_full
    riskdf$risk_Y_fraction <- (riskdf$risk_from_Y - full_risk)/no_info_Y
    riskdf$p_disagree_fraction <- riskdf$p_disagree/no_p_disagree
    
    
    return(riskdf)
}

SL_vimtest <- function(osl_fit, data, Xnodes, ...) {
    full_preds <- cv_predict_original(osl_fit)
    Y <- osl_fit$valY
    new_fit <- origami_SuperLearner(Y = Y, X = data[, Xnodes], folds = osl_fit$folds, family = osl_fit$fullFit$family, 
        SL.library = osl_fit$SL.library, method = osl_fit$fullFit$method, ...)
    new_preds <- cv_predict_original(new_fit)
    plot(full_preds, new_preds)
    browser()
    
}

#' @export
#'
tx_vim <- function(opt_obj, verbose = 2) {
    
    message_verbose("VIM for g", 1, verbose)
    obs <- with(opt_obj, SL_vim(fits$g_fit, data, nodes$Wnodes, cts.num = 5))
    obs$model <- "Observed"
    
    message_verbose("VIM for blip", 1, verbose)
    opt <- with(opt_obj, SL_vim(fits$blip_fit, data, nodes$Vnodes, split_preds = split_preds, cvfun = split_cv_SL, cts.num = 10, 
        blip_type = blip_type, maximize = maximize))
    opt$model <- "Optimal"
    
    vims <- rbind(obs, opt)
    
    return(vims)
} 
