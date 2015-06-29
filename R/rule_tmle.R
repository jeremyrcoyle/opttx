################################################################ TMLE for dynamic regimes

# TMLE for a rule works by redefining A=1 as 'followed rule' and A=0 as 'didn't follow rule' pA is updated accordingly
fit_ruletmle <- function(obsA, obsY, pA1, Q0W, Q1W, ruleA) {
    A <- as.numeric(obsA == ruleA)
    pA <- pA1 * ruleA + (1 - pA1) * (1 - ruleA)
    Qd <- Q1W * ruleA + Q0W * (1 - ruleA)
    tmledata <- data.frame(A = A, Y = obsY, gk = pA, Qk = Qd)
    gentmle(tmledata, ey1_estimate, ey1_update)
}

ruletmle <- function(obsA, obsY, pA1, Q0W, Q1W, ruleA) {
    results <- fit_ruletmle(obsA, obsY, pA1, Q0W, Q1W, ruleA)
    # put results in a dataframe for ease of usep
    est <- results$tmleests
    sd <- sqrt(results$ED2)/sqrt(length(obsA))
    lower <- est - qnorm(0.975) * sd
    upper <- est + qnorm(0.975) * sd
    data.frame(est = est, sd = sd, lower = lower, upper = upper)
}

# difference of two TMLEs using the delta method note that this may not work because sometimes the rules are nested
# which means in the limit one rule minus the other rule always has the same sign
two_ruletmle <- function(obsA, obsY, pA1, Q0W, Q1W, ruleA, rule2A) {
    results <- fit_ruletmle(obsA, obsY, pA1, Q0W, Q1W, ruleA)
    results2 <- fit_ruletmle(obsA, obsY, pA1, Q0W, Q1W, rule2A)
    # combined results and put in a dataframe for ease of usep
    est <- results$tmleests - results2$tmleests
    sd <- sqrt(mean((results$Dstar[[1]] - results2$Dstar[[1]])^2))/sqrt(length(obsA))
    lower <- est - qnorm(0.975) * sd
    upper <- est + qnorm(0.975) * sd
    data.frame(est = est, sd = sd, lower = lower, upper = upper)
}

# generates predictions from original dataset from origami_SuperLearner object
cv_predict_original <- function(osl_fit) {
    osl_fit$Z %*% osl_fit$coef
}

# extract the validation sets from split_preds for use in CV-TMLE
extract_val <- function(fold, split_preds) {
    v <- fold_index()
    valid_idx <- validation()
    
    val_preds <- sapply(split_preds, function(split_pred) {
        split_pred[[v]][valid_idx]
    })
    val_preds <- as.data.frame(val_preds)
    val_preds$index <- valid_idx
    result <- list(preds = val_preds)
    return(result)
}

all_ests <- function(data, folds, nodes, fits, split_preds, parallel = F) {
    # extract split-specific validation predictions
    val_preds <- cross_validate(extract_val, folds, split_preds, .parallel = parallel)$preds
    val_preds <- val_preds[order(val_preds$index), ]
    
    # fit static blip blipfit0=subopt(c('W2'),split_preds,folds,data,SL.library='SL.mean')
    cv_pblip <- cv_predict_original(fits$blip_fit)
    cv_optA <- as.numeric(cv_pblip > 0.5)
    A <- data[, nodes$Anode]
    Y <- data[, nodes$Ynode]
    ests <- with(val_preds, {
        est_0 <- ruletmle(A, Y, pA1, Q0W, Q1W, 0)
        est_1 <- ruletmle(A, Y, pA1, Q0W, Q1W, 1)
        est_A <- ruletmle(A, Y, pA1, Q0W, Q1W, A)
        est_opt <- ruletmle(A, Y, pA1, Q0W, Q1W, cv_optA)
        rbind(est_0, est_1, est_A, est_opt)
    })
    
    rules <- c("A=0", "A=1", "A=a", "A=d(V)")
    ests$rule <- rules
    
    # diffest <- tworuletmle2(data$A, data$Y, cv_pA1, cv_Q0W, cv_Q1W, old_optA, cv_optA)
    return(ests)
}


 
