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
    full_fit <- opt_obj$fits$rule_fit
    
    # fit rules on all subsets V_i\V
    newnodes <- opt_obj$nodes
    Vnodes <- newnodes$Vnodes
    system.time({
        reduced_fits <- lapply(Vnodes, function(drop_Vnode) {
            newnodes$Vnodes <- setdiff(opt_obj$nodes$Vnodes, drop_Vnode)
            reduced_fit <- with(opt_obj, learn_rule(data, folds, newnodes, split_preds, 
                val_preds, parallel = F, SL.library = SL.library, verbose = 3))
        })
    })
    testsamp <- data[, Vnodes]
    testsamp[, -1] <- testsamp[1, -1]
    
    system.time({
        stupid <- merge(data[, 1, drop = F], data[, -1])
        z <- get_test_preds(full_fit, stupid)
    })
    z
    get_cv_preds <- function(fit) {
        predict(fit, newdata = "cv-original", pred_fit = "joint", return_assignment = F)
    }
    preds_to_dV <- max.col
    EYd <- function(dV) {
        with(opt_obj, rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, 
            val_preds$QaW, dV))
    }
    EYd_comp <- function(reduced_dV, full_dV) {
        with(opt_obj, two_rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, 
            val_preds$QaW, full_dV, reduced_dV))
    }
    
    full_preds <- get_cv_preds(full_fit)
    full_dV <- preds_to_dV(full_preds)
    full_EYd <- EYd(dV)
    
    reduced_preds <- lapply(reduced_fits, get_cv_preds)
    reduced_dV <- lapply(reduced_preds, preds_to_dV)
    reduced_EYd_comp <- ldply(reduced_dV, EYd_comp, full_dV)
    n <- nrow(opt_obj$data)
    reduced_loglik <- -1 * sapply(reduced_preds, origami:::mn_loglik, full_dV, rep(1, 
        n))
    risk_01 <- function(x, y) {
        mean(x != y)
    }
    reduced_01 <- sapply(reduced_dV, risk_01, full_dV)
    reduced_dV_diff <- lapply(reduced_dV, `-`, full_dV)
    get_test_preds <- function(fit, test_data) {
        predict(fit, newdata = test_data[, fit$nodes$Vnodes], pred_fit = "joint", 
            return_assignment = F)
    }
    test_EYd <- function(dV, testdata) {
        mean(Qbar0(dV, testdata[, nodes$Wnodes]))
    }
    full_test_preds <- get_test_preds(full_fit, testdata)
    full_test_dV <- preds_to_dV(full_test_preds)
    full_test_EYd <- test_EYd(full_test_dV, testdata)
    
    
    reduced_test_preds <- lapply(reduced_fits, get_test_preds, testdata)
    reduced_test_dV <- lapply(reduced_test_preds, preds_to_dV)
    reduced_test_EYd <- sapply(reduced_test_dV, test_EYd, testdata)
    reduced_test_EYd_comp <- full_test_EYd - reduced_test_EYd
    reduced_EYd_comp$test_est <- reduced_test_EYd_comp
    reduced_test_01 <- sapply(reduced_test_dV, risk_01, full_test_dV)
    reduced_test_loglik <- -1 * sapply(reduced_test_preds, origami:::mn_loglik, full_test_dV, 
        rep(1, nrow(testdata)))
    list(preds = preds, dV = dV, EYd = EYd)
    diff_preds <- full$preds - reduced$preds
    diff_dV <- full$dV - reduced$dV
    diff_EYd <- with(opt_obj, two_rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], 
        val_preds$pA, val_preds$QaW, full$dV, reduced$dV))
    
    ggplot(long, aes(x = val, y = value, color = variable)) + facet_wrap(~node) + 
        geom_smooth(method = "loess") + geom_rug(alpha = 0.2, sides = "b")
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
