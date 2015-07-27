################################################################ TMLE for dynamic regimes

# TMLE for a rule works by redefining A=1 as 'followed rule' and A=0 as 'didn't
# follow rule' pA is updated accordingly todo: move to gentmle
fit_rule_tmle <- function(obsA, obsY, pA, QaW, dV) {
    Ad <- as.numeric(obsA == dV)
    dV_ind <- factor_to_indicators(dV)
    pAd <- rowSums(pA * dV_ind)
    Qd <- rowSums(QaW * dV_ind)
    tmledata <- data.frame(A = Ad, Y = obsY, gk = pAd, Qk = Qd)
    gentmle(tmledata, ey1_estimate, ey1_update)
}

rule_dripcw <- function(DR, dV) {
    dV_ind <- factor_to_indicators(dV)
    DR_dV <- rowSums(DR * dV_ind)
    est <- mean(DR_dV)
    sd <- sd(DR_dV)/sqrt(length(DR_dV))
    lower <- est - qnorm(0.975) * sd
    upper <- est + qnorm(0.975) * sd
    data.frame(est = est, sd = sd, lower = lower, upper = upper)
}

rule_tmle <- function(obsA, obsY, pA, QaW, dV) {
    results <- fit_rule_tmle(obsA, obsY, pA, QaW, dV)
    ci <- ci_gentmle(results)
    ci$parameter <- NULL
    
    return(ci)
}

# difference of two TMLEs using the delta method note that this may not work
# because sometimes the rules are nested which means in the limit one rule minus
# the other rule always has the same sign
two_rule_tmle <- function(obsA, obsY, pA, QaW, dV, dV2) {
    results <- fit_rule_tmle(obsA, obsY, pA, QaW, dV)
    results2 <- fit_rule_tmle(obsA, obsY, pA, QaW, dV2)
    # combined results and put in a dataframe for ease of use
    est <- results$tmleests - results2$tmleests
    sd <- sqrt(mean((results$Dstar[[1]] - results2$Dstar[[1]])^2))/sqrt(length(obsA))
    lower <- est - qnorm(0.975) * sd
    upper <- est + qnorm(0.975) * sd
    data.frame(est = est, sd = sd, lower = lower, upper = upper)
}


estimate_performance <- function(data, nodes, predictions, dV, perf_tmle = TRUE, 
    perf_dripcw = FALSE) {
    A <- data[, nodes$Anode]
    Y <- data[, nodes$Ynode]
    A_vals <- vals_from_factor(A)
    tmle_ests <- NULL
    
    if (perf_tmle) {
        static_ests <- ldply(A_vals, function(A_val) {
            rule_tmle(A, Y, predictions$pA, predictions$QaW, rep(A_val, length(Y)))
        })
        
        est_A <- rule_tmle(A, Y, predictions$pA, predictions$QaW, A)
        est_opt <- rule_tmle(A, Y, predictions$pA, predictions$QaW, dV)
        
        ests <- rbind(static_ests, est_A, est_opt)
        rules <- c(sprintf("A=%s", A_vals), "A=a_obs", "A=d(V)")
        ests$rule <- rules
        ests$estimator <- "TMLE"
        tmle_ests <- ests
    }
    
    dripcw_ests <- NULL
    if (perf_dripcw) {
        static_ests <- ldply(A_vals, function(A_val) {
            rule_dripcw(predictions$DR, rep(A_val, length(Y)))
        })
        
        est_A <- rule_dripcw(predictions$DR, A)
        est_opt <- rule_dripcw(predictions$DR, dV)
        
        ests <- rbind(static_ests, est_A, est_opt)
        rules <- c(sprintf("A=%s", A_vals), "A=a_obs", "A=d(V)")
        ests$rule <- rules
        ests$estimator <- "DR-IPCW"
        dripcw_ests <- ests
    }
    
    all_ests <- rbind(tmle_ests, dripcw_ests)
    
    return(all_ests)
} 
