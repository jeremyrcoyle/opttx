risk_from_osl <- function(osl_fit) {
    function(pred, Y) {
        method <- osl_fit$fullFit$method
        coef <- method$computeCoef(pred, Y, libraryNames = "split_SL", obsWeights = osl_fit$valWeights, 
            verbose = F)
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
    full_coefs <- osl_fit$coef
    which_nz <- which(full_coefs != 0)
    nz_coefs <- full_coefs[which_nz]
    nz_library <- osl_fit$SL.library[which_nz]
    
    riskdf <- ldply(seq_along(Xnodes), .progress = "text", function(idx, ...) {
        new_Xnodes <- Xnodes[-1 * idx]
        new_fit <- origami_SuperLearner(Y = as.vector(full_preds), X = data[, new_Xnodes], 
            folds = osl_fit$folds, family = "gaussian", SL.library = nz_library, 
            method = method.probloglik)
        
        new_preds <- new_fit$Z %*% nz_coefs  #cv_predict_original(new_fit)
        risk_from_full <- risk_fun(new_preds, full_preds)
        p_disagree <- mean((full_preds > 0.5) != (new_preds > 0.5))
        data.frame(node = Xnodes[idx], p_disagree = p_disagree, risk_from_full = risk_from_full)
    })
    
    no_fit <- origami_SuperLearner(Y = as.vector(full_preds), X = data[, Xnodes], 
        folds = osl_fit$folds, family = osl_fit$fullFit$family, SL.library = "SL.mean", 
        ...)
    no_pred <- cv_predict_original(no_fit)
    no_p_disagree <- mean((full_preds > 0.5) != (no_pred > 0.5))
    no_info_full <- risk_fun(no_pred, full_preds)
    no_info <- data.frame(node = "empty", p_disagree = no_p_disagree, risk_from_full = no_info_full)
    riskdf <- rbind(riskdf, no_info)
    riskdf$risk_fraction <- (riskdf$risk_from_full)/no_info_full
    riskdf$p_disagree_fraction <- riskdf$p_disagree/no_p_disagree
    
    
    return(riskdf)
}

SL_vimtest <- function(osl_fit, data, Xnodes, ...) {
    full_preds <- cv_predict_original(osl_fit)
    Y <- osl_fit$valY
    new_fit <- origami_SuperLearner(Y = Y, X = data[, Xnodes], folds = osl_fit$folds, 
        family = osl_fit$fullFit$family, SL.library = osl_fit$SL.library, method = osl_fit$fullFit$method, 
        ...)
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
    
    message_verbose("VIM for rule", 1, verbose)
    opt <- with(opt_obj, SL_vim(fits$blip_fit, data, nodes$Vnodes, split_preds = split_preds, 
        cvfun = split_cv_SL, cts.num = 10, blip_type = blip_type, maximize = maximize))
    opt$model <- "Optimal"
    
    vims <- rbind(obs, opt)
    
    return(vims)
}

sub_rule <- function(opt_obj, V) {
    opt_obj <- result
    rule_fit <- opt_obj$fits$rule_fit
    val_preds <- with(opt_obj, extract_vals(folds, split_preds))
    
    # invert -- calc all parameters for each reduced, then extract from lists write
    # calc parameter function so can also apply to full
    newnodes <- opt_obj$nodes
    Vnodes <- newnodes$Vnodes
    reduced_fits <- lapply(Vnodes, function(drop_Vnode) {
        newnodes$Vnodes <- setdiff(opt_obj$nodes$Vnodes, drop_Vnode)
        reduced_fit <- with(opt_obj, learn_rule(data, folds, newnodes, split_preds, 
            val_preds, parallel = F, SL.library = SL.library, verbose = 3))
    })
    
    full_preds <- predict(rule_fit, newdata = "cv-original", pred_fit = "joint", 
        return_assignment = F)
    reduced_preds <- lapply(reduced_fits, predict, newdata = "cv-original", pred_fit = "joint", 
        return_assignment = F)
    diff_dfs <- lapply(reduced_preds, function(reduced_pred) {
        data.frame(diff = reduced_pred - full_preds)
    })
    diff_dfs <- mapply(function(df, node) {
        df$node <- node
        df$val <- data[, node]
        df
    }, df = diff_dfs, node = Vnodes, SIMPLIFY = F)
    diff_df <- do.call(rbind, diff_dfs)
    head(diff_df)
    long <- melt(diff_df, id = c("node", "val"))
    
    head(long)
    ggplot(long, aes(x = val, y = value, color = variable)) + facet_wrap(~node) + 
        geom_smooth(method = "loess") + geom_rug(alpha = 0.2, sides = "b")
    full_dV <- max.col(full_preds)
    reduced_dV <- lapply(reduced_preds, max.col)
    
    full_EYd <- with(opt_obj, rule_tmle(data$A, data$Y, val_preds$pA, val_preds$QaW, 
        full_dV))
    reduced_EYd <- with(opt_obj, ldply(reduced_dV, function(dV) two_rule_tmle(data$A, 
        data$Y, val_preds$pA, val_preds$QaW, full_dV, dV)))
    reduced_EYd
    test_reduced_dV <- lapply(reduced_fits, function(fit) predict(fit, newdata = testdata[, 
        fit$nodes$Vnodes], pred_fit = "joint"))
    reduced_E0Yd <- colMeans(sapply(test_reduced_dV, Qbar0, testdata[, newnodes$Wnodes]))
    test_full_dV <- predict(rule_fit, newdata = testdata[, Vnodes], pred_fit = "joint")
    z <- sapply(reduced_dV, function(pred_dV) mean(full_dV == pred_dV))
    z0 <- sapply(test_reduced_dV, function(pred_dV) mean(test_full_dV == pred_dV))
    plot(z0, z)
    full_E0Yd <- mean(Qbar0(test_full_dV, testdata[, newnodes$Wnodes]))
    diff0 <- full_E0Yd - reduced_E0Yd
    plot(diff0, reduced_EYd$est)
    newnodes$Vnodes <- setdiff(opt_obj$nodes$Vnodes, "W5")
    rule_fit3 <- with(opt_obj, learn_rule(data, folds, newnodes, split_preds, val_preds, 
        parallel = F, SL.library = SL.library, verbose = 3))
    
    cv_pdV <- with(rule_fit, predict(joint_fit, QaV_fits, class_fit, newdata = "cv-original"))
    
    cv_pdV3 <- with(rule_fit3, predict(joint_fit, QaV_fits, class_fit, newdata = "cv-original"))
    quantile(cv_pdV - cv_pdV2)
    quantile(cv_pdV - cv_pdV3)
    
    
    library(reshape2)
    
    ggplot(dropdata, aes(x = var, y = value, color = variable)) + geom_smooth(method = "loess") + 
        facet_wrap(~name) + theme_bw()
}
test <- function(opt_obj) {
    newnodes <- opt_obj$nodes
    Vnodes <- newnodes$Vnodes
    drop_Vnode <- Vnodes[1]
    cv_dV <- predict(rule_fit, newdata = "cv-original", pred_fit = "joint")
    full_fit <- multinomial_SuperLearner(cv_dV, data[, Vnodes], folds = opt_obj$folds)
    full_pred <- predict(full_fit, newdata = "cv-original")$pred
    
    dropdata2 <- ldply(Vnodes, function(drop_Vnode) {
        newnodes$Vnodes <- setdiff(opt_obj$nodes$Vnodes, drop_Vnode)
        red_fit <- multinomial_SuperLearner(cv_dV, data[, newnodes$Vnodes], folds = opt_obj$folds)
        red_pred <- predict(red_fit, newdata = "cv-original")$pred
        diffdf <- data.frame(diff = full_pred - red_pred, var = data[, drop_Vnode], 
            name = drop_Vnode)
        long <- melt(diffdf, id = c("var", "name"))
    })
    
    
    
    ggplot(dropdata2, aes(x = var, y = value, color = variable)) + geom_point(alpha = 0.1) + 
        geom_smooth(method = "loess") + facet_wrap(~name) + theme_bw()
    test2 <- multinomial_SuperLearner(cv_dV, data[, Wnodes], folds = opt_obj$folds)
}

test <- function(opt_obj) {
    
    newnodes <- opt_obj$nodes
    Vnodes <- newnodes$Vnodes
    drop_Vnode <- Vnodes[1]
    cv_dV <- predict(rule_fit, newdata = "cv-original", pred_fit = "joint")
    full_fit <- multinomial_SuperLearner(cv_dV, data[, Vnodes], folds = opt_obj$folds)
    full_fit <- opt_obj$fits$g_fit
    
    outcome <- full_fit$valY
    
    full_pred <- predict(full_fit, newdata = "cv-original")$pred
    full_cat <- max.col(full_pred)
    full_loglik <- mn_loglik(full_pred, factor_to_indicators(outcome), full_fit$valWeights)
    
    full_acc <- mean(full_cat == outcome)
    reduced_fits <- lapply(Vnodes, function(drop_Vnode) {
        newnodes$Vnodes <- setdiff(opt_obj$nodes$Vnodes, drop_Vnode)
        red_fit <- multinomial_SuperLearner(outcome, data[, newnodes$Vnodes], folds = opt_obj$folds)
    })
    
    
    reduced_preds <- lapply(reduced_fits, function(fit) predict(fit, newdata = "cv-original")$pred)
    reduced_cats <- sapply(reduced_preds, max.col)
    reduced_acc <- apply(reduced_cats, 2, function(x) mean(x == outcome))
    (full_acc - reduced_acc)/reduced_acc
    
    reduced_logliks <- sapply(reduced_preds, mn_loglik, factor_to_indicators(outcome), 
        full_fit$valWeights)
    lr <- 2 * (full_loglik - reduced_logliks)
    pchisq(lr, 10)
    red_fit$cvRisk
    ggplot(dropdata3, aes(x = var, y = value, color = variable)) + geom_point(alpha = 0.1) + 
        geom_smooth(method = "loess") + facet_wrap(~name) + theme_bw()
    test2 <- multinomial_SuperLearner(cv_dV, data[, Wnodes], folds = opt_obj$folds)
} 
