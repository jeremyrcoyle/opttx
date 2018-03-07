
wald <- function(x) {
    est <- mean(x)
    se <- sd(x)/sqrt(length(x))
    lower <- est - qnorm(0.975) * se
    upper <- est + qnorm(0.975) * se
    c(est = est, se = se, lower = lower, upper = upper)
}

get_cv_preds <- function(fit) {
    predict(fit, newdata = "cv-original")$pred
}

EYd <- function(fit_dV, opt_obj) {
    with(opt_obj, rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, val_preds$QaW, 
        fit_dV))
}

EYd_comp <- function(reduced_dV, full_dV, opt_obj) {
    with(opt_obj, two_rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, val_preds$QaW, 
        full_dV, reduced_dV))
}

risk_mse <- function(x, y) {
    loss <- (x - y)^2
    avg_loss <- apply(loss, 1, mean)
    wald(avg_loss)
}

risk_diff_mse <- function(reduced_x, full_x, y) {
    loss <- (reduced_x - y)^2 - (full_x - y)^2
    avg_loss <- apply(loss, 1, mean)
    wald(avg_loss)
}

risk_01 <- function(x, y) {
    wald(x != y)
}

risk_diff_01 <- function(reduced_x, full_x, y) {
    wald((reduced_x != y) - (full_x != y))
}

get_test_preds <- function(fit, testdata) {
    predict(fit, newdata = testdata[, fit$Vnodes])$pred
}

test_EYd <- function(dV, testdata, Qbar0) {
    mean(Qbar0(dV, testdata))
}

backward_vim <- function(opt_obj, testdata = NULL, Qbar0 = NULL) {
    # opt_obj <- result
    full_fit <- opt_obj$fits$rule_fit
    
    # fit rules on all subsets V_i\V
    Vnodes <- opt_obj$nodes$Vnodes
    full_fit$Vnodes <- Vnodes
    rule_args <- opt_obj$rule_args
    
    full_fit <- split_to_nested(full_fit, rule_args)
    
    reduced_fits <- lapply(Vnodes, function(drop_Vnode) {
        cat(sprintf("Vnode %s\n", drop_Vnode))
        newV <- setdiff(Vnodes, drop_Vnode)
        rule_args$X <- opt_obj$data[, newV, drop = F]
        reduced_fit <- split_from_args(rule_args)
        reduced_fit <- split_to_nested(reduced_fit, rule_args)
        reduced_fit$Vnodes <- newV
        
        reduced_fit
    })
    
    A_vals <- vals_from_factor(opt_obj$val_preds$A)
    
    
    full_preds <- get_cv_preds(full_fit)
    reduced_preds <- lapply(reduced_fits, get_cv_preds)
    # wait this is the DR-IPCW mapping for EY_a
    true_preds <- full_fit$valY
    
    reduced_mse <- ldply(reduced_preds, risk_diff_mse, full_preds, true_preds)
    reduced_mse$Vnode <- Vnodes
    reduced_mse$metric <- "regression mse"
    
    
    full_dV <- dV_from_preds(full_preds, A_vals)
    reduced_dV <- lapply(reduced_preds, dV_from_preds, A_vals)
    true_dV <- dV_from_preds(true_preds, A_vals)
    
    ldply(reduced_dV, risk_01, true_dV)
    risk_01(full_dV, true_dV)
    reduced_01 <- ldply(reduced_dV, risk_diff_01, full_dV, true_dV)
    reduced_01$Vnode <- Vnodes
    reduced_01$metric <- "rule disagreement"
    
    # full_EYd <- EYd(full_dV, opt_obj) reduced_EYd <- ldply(reduced_dV, EYd, opt_obj)
    reduced_EYd_comp <- ldply(reduced_dV, EYd_comp, full_dV, opt_obj)
    reduced_EYd_comp$Vnode <- Vnodes
    reduced_EYd_comp$metric <- "value difference"
    
    vimdf <- rbind.fill(reduced_mse, reduced_01, reduced_EYd_comp)
    
    Vnode_order <- reduced_mse$Vnode[order(reduced_mse$est)]
    vimdf$Vnode <- factor(vimdf$Vnode, levels = Vnode_order)
    
    # pdf('vim.pdf', height = 4, width = 6) ggplot(vimdf, aes(y = Vnode, x = est, xmin = lower, xmax =
    # upper)) + geom_point() + geom_errorbarh() + facet_wrap(~metric, scales = 'free') + theme_bw()
    # dev.off()
    
    
    if (!is.null(testdata)) {
        full_test_preds <- get_test_preds(full_fit, testdata)
        reduced_test_preds <- lapply(reduced_fits, get_test_preds, testdata)
        # shoudl be testdata Y
        true_test_preds <- testdata$Q0aW
        
        reduced_test_mse <- ldply(reduced_test_preds, risk_diff_mse, full_test_preds, true_test_preds)
        test_mse_df <- data.frame(test = reduced_test_mse$est, Vnode = Vnodes, metric = "regression mse")
        
        full_test_dV <- dV_from_preds(full_test_preds, A_vals)
        reduced_test_dV <- lapply(reduced_test_preds, dV_from_preds, A_vals)
        true_test_dV <- dV_from_preds(true_test_preds, A_vals)
        
        
        reduced_test_01 <- ldply(reduced_test_dV, risk_diff_01, full_test_dV, true_test_dV)
        test_01_df <- data.frame(test = reduced_test_01$est, Vnode = Vnodes, metric = "rule disagreement")
        
        full_test_EYd <- test_EYd(full_test_dV, testdata, Qbar0)
        reduced_test_EYd <- sapply(reduced_test_dV, test_EYd, testdata, Qbar0)
        
        reduced_test_EYd_comp <- full_test_EYd - reduced_test_EYd
        
        test_EYd_comp_df <- data.frame(test = reduced_test_EYd_comp, Vnode = Vnodes, metric = "value difference")
        
        testdf <- rbind(test_mse_df, test_01_df, test_EYd_comp_df)
        
        vimdf <- merge(vimdf, testdf)
        vimdf$covered <- vimdf$lower < vimdf$test & vimdf$test < vimdf$upper
    }
    
    vimresult <- list(full_preds = full_preds, reduced_preds = reduced_preds, vimdf = vimdf)
    
    return(vimresult)
}


backward_vim2 <- function(opt_obj, testdata = NULL, Qbar0 = NULL) {
    # opt_obj <- result
    full_fit <- opt_obj$fits$rule_fit
    
    # fit rules on all subsets V_i\V
    Vnodes <- opt_obj$nodes$Vnodes
    
    rule_args <- opt_obj$rule_args
    
    full_fit <- split_from_args(rule_args)
    full_fit$Vnodes <- Vnodes
    reduced_fits <- lapply(Vnodes, function(drop_Vnode) {
        cat(sprintf("Vnode %s\n", drop_Vnode))
        newV <- setdiff(Vnodes, drop_Vnode)
        rule_args$X <- opt_obj$data[, newV, drop = F]
        reduced_fit <- split_from_args(rule_args)
        # reduced_fit <- split_to_nested(reduced_fit, rule_args)
        reduced_fit$Vnodes <- newV
        
        reduced_fit
    })
    
    A_vals <- vals_from_factor(opt_obj$val_preds$A)
    
    
    full_preds <- get_cv_preds(full_fit)
    reduced_preds <- lapply(reduced_fits, get_cv_preds)
    true_preds <- full_fit$valY
    
    reduced_mse <- ldply(reduced_preds, risk_diff_mse, full_preds, true_preds)
    reduced_mse$Vnode <- Vnodes
    reduced_mse$metric <- "regression mse"
    
    full_dV <- dV_from_preds(full_preds, A_vals)
    reduced_dV <- lapply(reduced_preds, dV_from_preds, A_vals)
    true_dV <- dV_from_preds(true_preds, A_vals)
    
    ldply(reduced_dV, risk_01, true_dV)
    risk_01(full_dV, true_dV)
    reduced_01 <- ldply(reduced_dV, risk_diff_01, full_dV, true_dV)
    reduced_01$Vnode <- Vnodes
    reduced_01$metric <- "rule disagreement"
    
    # full_EYd <- EYd(full_dV, opt_obj) reduced_EYd <- ldply(reduced_dV, EYd, opt_obj)
    reduced_EYd_comp <- ldply(reduced_dV, EYd_comp, full_dV, opt_obj)
    reduced_EYd_comp$Vnode <- Vnodes
    reduced_EYd_comp$metric <- "value difference"
    
    vimdf <- rbind.fill(reduced_mse, reduced_01, reduced_EYd_comp)
    
    Vnode_order <- reduced_mse$Vnode[order(reduced_mse$est)]
    vimdf$Vnode <- factor(vimdf$Vnode, levels = Vnode_order)
    
    # pdf('vim.pdf', height = 4, width = 6) ggplot(vimdf, aes(y = Vnode, x = est, xmin = lower, xmax =
    # upper)) + geom_point() + geom_errorbarh() + facet_wrap(~metric, scales = 'free') + theme_bw()
    # dev.off()
    
    
    if (!is.null(testdata)) {
        full_test_preds <- get_test_preds(full_fit, testdata)
        reduced_test_preds <- lapply(reduced_fits, get_test_preds, testdata)
        true_test_preds <- testdata$Q0aW
        
        reduced_test_mse <- ldply(reduced_test_preds, risk_diff_mse, full_test_preds, true_test_preds)
        test_mse_df <- data.frame(test = reduced_test_mse$est, Vnode = Vnodes, metric = "regression mse")
        
        full_test_dV <- dV_from_preds(full_test_preds, A_vals)
        reduced_test_dV <- lapply(reduced_test_preds, dV_from_preds, A_vals)
        true_test_dV <- dV_from_preds(true_test_preds, A_vals)
        
        
        reduced_test_01 <- ldply(reduced_test_dV, risk_diff_01, full_test_dV, true_test_dV)
        test_01_df <- data.frame(test = reduced_test_01$est, Vnode = Vnodes, metric = "rule disagreement")
        
        full_test_EYd <- test_EYd(full_test_dV, testdata, Qbar0)
        reduced_test_EYd <- sapply(reduced_test_dV, test_EYd, testdata, Qbar0)
        
        reduced_test_EYd_comp <- full_test_EYd - reduced_test_EYd
        
        test_EYd_comp_df <- data.frame(test = reduced_test_EYd_comp, Vnode = Vnodes, metric = "value difference")
        
        testdf <- rbind(test_mse_df, test_01_df, test_EYd_comp_df)
        
        vimdf <- merge(vimdf, testdf)
        vimdf$covered <- vimdf$lower < vimdf$test & vimdf$test < vimdf$upper
    }
    
    vimresult <- list(full_preds = full_preds, reduced_preds = reduced_preds, vimdf = vimdf)
    
    return(vimresult)
}

vis_preds <- function(vimresult, data) {
    index <- 8
    SL.library <- c("SL.glm", "SL.loess", "SL.mean")
    mvSL.library <- sl_to_mv_library(SL.library)
    all_diffs <- ldply(setdiff(1:length(Vnodes), 7), function(index) {
        
        diffmat <- vimresult$reduced_preds[[index]]  #vimresult$full_preds - vimresult$reduced_preds[[index]]
        Vnode <- Vnodes[index]
        print(Vnode)
        V <- data[, Vnode, drop = F]
        
        # sl <- origami_SuperLearner(folds = result$folds, Y = diffmat, X = V, SL.library = mvSL.library,
        # method = method.mvSL(method.NNLS())) predmat <- predict(sl, newdata = V)$pred
        diffs <- as.data.frame(diffmat)
        diffs$Vnode <- Vnode
        diffs$V <- V[, 1]
        diffs
        
    })
    
    A_vals <- vals_from_factor(data[, Anode])
    
    names(all_diffs)[1:3] <- as.character(A_vals)
    diffdf <- melt(all_diffs, id = c("Vnode", "V"))
    diffdf
    
    # testdiff_preds <- compare_preds(full_test_preds, reduced_test_preds, testdata)
    ct.num <- sapply(data[, Vnodes], function(x) length(unique(x)))
    # diff_preds <- diff_preds[(diff_preds$Vnode %in% goodnodes), ]
    V_facts <- Vnodes[which(ct.num < 5)]
    diff_conts <- diffdf[!(diffdf$Vnode %in% V_facts), ]
    
    # trim extreme values so we show only were there's support
    diff_conts <- ddply(diff_conts, .(Vnode), function(nodedata) {
        quants <- quantile(nodedata$V, c(0.025, 0.975))
        nodedata[quants[1] < nodedata$V & nodedata$V < quants[2], ]
    })
    
    pdf("vim_cont.pdf", height = 6, width = 8)
    ggplot(diff_conts, aes(x = V)) + facet_wrap(~Vnode, scales = "free") + geom_point(aes(y = value, 
        color = variable, group = variable), alpha = 0.2) + geom_smooth(aes(y = value, color = variable), 
        method = "loess", se = F) + geom_rug(data = diff_conts[diff_conts$variable == diff_conts$variable[1], 
        ], color = "black", alpha = 0.2, sides = "b") + theme_bw()
    dev.off()
    
    pdf("~/Dropbox/thesis/quals/slides/graphs/vim_pulse.pdf", height = 3, width = 6)
    ggplot(diff_conts[diff_conts$Vnode == "PULSEVS", ], aes(x = V)) + geom_line(aes(y = value, color = variable, 
        group = variable)) + geom_rug(data = diff_conts[diff_conts$variable == diff_conts$variable[1] & 
        diff_conts$Vnode == "PULSEVS", ], color = "black", alpha = 0.2, sides = "b") + theme_bw() + 
        xlab("Pulse (BPM)") + ylab("Q(V)-Q(V')")
    dev.off()
    
    dodge <- position_dodge(width = 0.5)
    injtype <- diffdf[diffdf$Vnode == "X_1STPENE", ]
    injtype$V <- factor(injtype$V, levels = c(1, 2), labels = c("Blunt", "Penetrating"))
    pdf("~/Dropbox/thesis/quals/slides/graphs/vim_injurtytype.pdf", height = 3, width = 4)
    ggplot(injtype, aes(x = factor(V))) + geom_point(aes(y = value, color = variable), position = dodge) + 
        theme_bw() + ylab("Q(V)-Q(V')") + xlab("Injury Type")
    dev.off()
    aggregate(value ~ Vnode, diff_preds, var)
    aggregate(value ~ Vnode, testdiff_preds, var)
    diff_dV <- full$dV - reduced$dV
    diff_EYd <- with(opt_obj, two_rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, 
        val_preds$QaW, full$dV, reduced$dV))
}

vis_preds <- function(vimresult, data) {
    index <- 8
    mnSL.library <- c("mvSL.multinom", "mnSL.mean")
    
    all_diffs <- ldply(1:length(Vnodes), function(index) {
        full_dV <- dV_from_preds(vimresult$full_preds, A_vals)
        red_dV <- dV_from_preds(vimresult$reduced_preds[[index]])
        Vnode <- Vnodes[index]
        V <- data[, Vnode, drop = F]
        
        full_sl <- origami_SuperLearner(folds = result$folds, Y = full_dV, X = V, SL.library = mnSL.library, 
            method = method.mnNNloglik(), family = list(family = "multinomial"))
        red_sl <- origami_SuperLearner(folds = result$folds, Y = red_dV, X = V, SL.library = mnSL.library, 
            method = method.mnNNloglik(), family = list(family = "multinomial"))
        full_pred <- predict(full_sl, newdata = V)$pred
        red_pred <- predict(red_sl, newdata = V)$pred
        
        diffs <- as.data.frame(full_pred - red_pred)
        diffs$Vnode <- Vnode
        diffs$V <- V[, 1]
        diffs
        
    })
    
    A_vals <- vals_from_factor(data[, Anode])
    
    names(all_diffs)[1:3] <- as.character(A_vals)
    diffdf <- melt(all_diffs, id = c("Vnode", "V"))
    
    
    # testdiff_preds <- compare_preds(full_test_preds, reduced_test_preds, testdata)
    ct.num <- sapply(data[, Vnodes], function(x) length(unique(x)))
    # diff_preds <- diff_preds[(diff_preds$Vnode %in% goodnodes), ]
    V_facts <- Vnodes[which(ct.num < 5)]
    diff_conts <- diffdf[!(diffdf$Vnode %in% V_facts), ]
    
    # trim extreme values so we show only were there's support
    diff_conts <- ddply(diff_conts, .(Vnode), function(nodedata) {
        quants <- quantile(nodedata$V, c(0.025, 0.975))
        nodedata[quants[1] < nodedata$V & nodedata$V < quants[2], ]
    })
    
    pdf("vim_cont.pdf", height = 6, width = 8)
    ggplot(diff_conts, aes(x = V)) + facet_wrap(~Vnode, scales = "free") + geom_line(aes(y = value, 
        color = variable, group = variable)) + geom_rug(data = diff_conts[diff_conts$variable == diff_conts$variable[1], 
        ], color = "black", alpha = 0.2, sides = "b") + theme_bw()
    dev.off()
    
    pdf("~/Dropbox/thesis/quals/slides/graphs/vim_pulse.pdf", height = 3, width = 6)
    ggplot(diff_conts[diff_conts$Vnode == "PULSEVS", ], aes(x = V)) + geom_line(aes(y = value, color = variable, 
        group = variable)) + geom_rug(data = diff_conts[diff_conts$variable == diff_conts$variable[1] & 
        diff_conts$Vnode == "PULSEVS", ], color = "black", alpha = 0.2, sides = "b") + theme_bw() + 
        xlab("Pulse (BPM)") + ylab("P(d(V)=a)-P(d(V')=a)")
    dev.off()
    
    dodge <- position_dodge(width = 0.5)
    injtype <- diffdf[diffdf$Vnode == "X_1STPENE", ]
    injtype$V <- factor(injtype$V, levels = c(1, 2), labels = c("Blunt", "Penetrating"))
    pdf("~/Dropbox/thesis/quals/slides/graphs/vim_injurtytype.pdf", height = 3, width = 4)
    ggplot(injtype, aes(x = factor(V))) + geom_point(aes(y = value, color = variable), position = dodge) + 
        theme_bw() + ylab("Q(V)-Q(V')") + xlab("Injury Type")
    dev.off()
    aggregate(value ~ Vnode, diff_preds, var)
    aggregate(value ~ Vnode, testdiff_preds, var)
    diff_dV <- full$dV - reduced$dV
    diff_EYd <- with(opt_obj, two_rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, 
        val_preds$QaW, full$dV, reduced$dV))
}


# todo: be clever about categorical V
forward_vim <- function(opt_obj, V) {
    opt_obj <- result
    full_fit <- opt_obj$fits$rule_fit
    
    # fit rules on all subsets V_i\V
    Vnodes <- opt_obj$nodes$Vnodes
    full_fit$Vnodes <- Vnodes
    rule_args <- opt_obj$rule_args
    A_vals <- vals_from_factor(opt_obj$val_preds$A)
    
    EYd <- function(fit_dV) {
        with(opt_obj, rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, val_preds$QaW, 
            fit_dV))
    }
    EYd_comp <- function(reduced_dV, full_dV) {
        with(opt_obj, two_rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, val_preds$QaW, 
            full_dV, reduced_dV))
    }
    
    
    rule_args$Y <- full_preds
    rule_split_preds <- opt_obj$split_preds
    rule_split_preds$DR <- with(opt_obj, cross_validate(opttx_refit_split_preds, folds, data, nodes, 
        fits, .combine = F))$DR
    rule_full_preds <- opt_obj$full_preds
    rule_full_preds$DR <- with(opt_obj, opttx_refit_split_preds(folds[[1]], data, nodes, fits, use_full = T))$DR
    # rule_val_preds <- extract_vals(folds, rule_split_preds)
    
    rule_args$split_preds <- rule_split_preds
    rule_args$full_preds <- rule_full_preds
    rule_args$Y <- rule_full_preds$DR
    rule_args$SL.library <- setdiff(rule_args$SL.library, "mvSL.glmnet")  #drop glmnets
    # be smarter about how we combine rule_fit for different cats (right now it's a stupid average)
    drop_Vnode <- Vnodes[[1]]
    system.time({
        reduced_fits <- lapply(Vnodes, function(drop_Vnode) {
            
            cat(sprintf("Vnode %s\n", drop_Vnode))
            # newV <- setdiff(Vnodes, drop_Vnode)
            
            rule_args$X <- opt_obj$data[, drop_Vnode, drop = F]
            reduced_fit <- split_from_args(rule_args)
            reduced_fit$Vnodes <- drop_Vnode
            
            reduced_fit
        })
    })
    
    rule_args$X <- opt_obj$data[, Vnodes, drop = F]
    
    
    A_vals <- vals_from_factor(opt_obj$val_preds$A)
    null_preds <- predict(reduced_fits[[1]], newdata = "cv-original")$library[, , "mvSL.mean"]
    null_dV <- dV_from_preds(null_preds, A_vals)
    null_EYd <- EYd(null_dV)
    
    
    reduced_preds <- lapply(reduced_fits, get_cv_preds)
    reduced_dV <- lapply(reduced_preds, dV_from_preds, A_vals)
    reduced_EYd <- ldply(reduced_dV, EYd)
    reduced_EYd_comp <- ldply(reduced_dV, function(reduced) EYd_comp(null_dV, reduced))
    reduced_01 <- sapply(reduced_dV, risk_01, null_dV)
    reduced_preds_diff <- lapply(reduced_preds, `-`, null_preds)
    reduced_mse <- sapply(reduced_preds_diff, function(x) var(as.vector(x)))
    
    reduced_EYd_comp$Vnode <- Vnodes
    reduced_EYd_comp$mse <- reduced_mse
    reduced_EYd_comp$pdisagree <- reduced_01
    reduced <- reduced_EYd_comp[order(reduced_EYd_comp$mse), ]
    
    get_test_preds <- function(fit, testdata) {
        predict(fit, newdata = testdata[, fit$Vnodes])$pred
    }
    
    test_EYd <- function(dV, testdata) {
        mean(Qbar0(dV, testdata))
    }
    
    full_test_preds <- get_test_preds(full_fit, testdata)
    full_test_dV <- dV_from_preds(full_test_preds)
    full_test_EYd <- test_EYd(full_test_dV, testdata)
    
    
    reduced_test_preds <- lapply(reduced_fits, get_test_preds, testdata)
    reduced_test_dV <- lapply(reduced_test_preds, dV_from_preds)
    reduced_test_EYd <- sapply(reduced_test_dV, test_EYd, testdata)
    reduced_EYd$test_est <- reduced_test_EYd
    reduced_test_EYd_comp <- full_test_EYd - reduced_test_EYd
    reduced_EYd_comp$test_est <- reduced_test_EYd_comp
    reduced_test_01 <- sapply(reduced_test_dV, risk_01, full_test_dV)
    reduced_test_preds_diff <- lapply(reduced_test_preds, `-`, full_test_preds)
    reduced_test_mse <- sapply(reduced_test_preds_diff, function(x) var(as.vector(x)))
    
    plot(reduced_test_mse, reduced_mse)
    list(preds = preds, dV = dV, EYd = EYd)
    
    
    compare_preds <- function(full, reduced, Vdata) {
        all_diffs <- ldply(1:length(Vnodes), function(index) {
            # diffs <- as.data.frame(full - reduced[[index]])
            diffs <- as.data.frame(reduced[[index]])
            diffs$Vnode <- Vnodes[index]
            diffs$V <- Vdata[, Vnodes[index]]
            diffs
        })
        long <- melt(all_diffs, id = c("Vnode", "V"))
    }
    
    plot_reduced_preds <- lapply(reduced_fits, function(x) predict(x, newdata = opt_obj$data)$pred)
    plot_null_preds <- predict(reduced_fits[[1]], newdata = opt_obj$data)$library[, , 6]
    diff_preds <- compare_preds(plot_null_preds, plot_reduced_preds, opt_obj$data)
    # diff_preds <- compare_preds(factor_to_indicators(full_dV),
    # lapply(reduced_dV,factor_to_indicators), opt_obj$data)
    diff_preds$variable <- levels(opt_obj$val_preds$A)[diff_preds$variable]
    # testdiff_preds <- compare_preds(full_test_preds, reduced_test_preds, testdata)
    goodnodes <- tail(levels(reducedlong$Vnode), 100)
    ct.num <- sapply(opt_obj$data[, Vnodes], function(x) length(unique(x)))
    diff_preds <- diff_preds[(diff_preds$Vnode %in% goodnodes), ]
    V_facts <- Vnodes[which(ct.num < 10)]
    diff_conts <- diff_preds[!(diff_preds$Vnode %in% V_facts), ]
    test <- ddply(diff_conts, .(Vnode), function(nodedata) {
        quants <- quantile(nodedata$V, c(0.025, 0.975))
        print(nodedata$Vnode[1])
        print(quants)
        nodedata[quants[1] < nodedata$V & nodedata$V < quants[2], ]
    })
    
    dim(diff_conts)
    pdf("vim_cont.pdf", height = 6, width = 8)
    ggplot(test, aes(x = V)) + facet_wrap(~Vnode, scales = "free") + geom_line(aes(y = value, color = variable, 
        group = variable)) + geom_rug(data = test[test$variable == test$variable[1], ], color = "black", 
        alpha = 0.8, sides = "b") + theme_bw() + coord_cartesian(ylim = c(-0.1, 0.1))
    dev.off()
    diff_facts <- diff_preds[diff_preds$Vnode %in% V_facts, ]
    ggplot(diff_facts, aes(x = factor(V))) + facet_wrap(~Vnode, scales = "free") + geom_point(aes(y = value, 
        color = variable), position = "dodge") + theme_bw() + coord_cartesian(ylim = range(diff_facts$value))
    aggregate(value ~ Vnode, diff_preds, var)
    aggregate(value ~ Vnode, testdiff_preds, var)
    diff_dV <- full$dV - reduced$dV
    diff_EYd <- with(opt_obj, two_rule_tmle(data[, nodes$Anode], data[, nodes$Ynode], val_preds$pA, 
        val_preds$QaW, full$dV, reduced$dV))
    # do this for the estimates after TMLE!
    list(full_preds = full_preds, reudced_preds = reduced_preds, vimdf = vimdf)
    
} 
