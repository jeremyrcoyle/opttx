#' @title Example library specification
#' @description Separate libraries must be specified for Q, g, QaV, and the classifier
opt_tmle.SL.library <- list(Q = c("SL.glm", "SL.glmem", "SL.glmnetprob", "SL.step.forward", 
    "SL.gam", "SL.rpart", "SL.rpartPrune", "SL.mean"), g = c("mnSL.randomForest", 
    "mnSL.glmnet", "mnSL.multinom", "mnSL.mean"), QaV = c("SL.glm", "SL.glmnetprob", 
    "SL.step.forward", "SL.gam", "SL.rpart", "SL.rpartPrune", "SL.mean"), class = c("mnSL.randomForest", 
    "mnSL.glmnet", "mnSL.multinom", "mnSL.mean"))

#' @title opt_tmle
#' @description Estimation of the Optimal Treatment rule using Super Learner and mean performance using CV-TMLE
#' To avoid nesting cross-validation, it uses split-specfic estimates of Q and g to estimate the rule, and 
#' 'split-specific' estimates of the rule in CV-TMLE to estimate mean performance
#' @param data data.frame containing the relevant variable
#' @param Wnodes, vector of column names indicating covariates
#' @param Anode, column name of treatment
#' @param Ynode, column name of outcome
#' @param Vnodes, vector of column names to base the treatment on
#' @param stratifyAY, logical: should we stratify the cross-validation based on (A,Y) pairs
#' @param SL_library, list of SuperLearner libraries for the various models. See \code{\link{opt_tmle.SL_library}} for an example.
#' @param verbose, integer that controls the verbosity of the output (higher is more verbose)
#' @param parallel, logical: should foreach parallelization be used?
#' @example /inst/examples/opttx.R
#' 
#' @export
opt_tmle <- function(data, Wnodes = grep("^W", names(data), value = T), Anode = "A", 
    Ynode = "Y", Vnodes = Wnodes, stratifyAY = TRUE, SL.library = opt_tmle.SL.library, 
    verbose = 3, parallel = F, perf_tmle = T, perf_dripcw = F, perf_cv = T, perf_full = F, 
    ...) {
    
    # ensure A is a factor
    data[, Anode] <- as.factor(data[, Anode])
    
    if (stratifyAY) {
        AYstrata <- sprintf("%s %s", data[, Anode], data[, Ynode])
        folds <- make_folds(strata_ids = AYstrata, V = 10)
    } else {
        folds <- make_folds(V = 10)
    }
    
    # if (!maximize) { todo:reimplement inversion of Y for adverse outcomes!  perhaps
    # easiest in opttx_split_preds }
    
    # possibly we should make these lists the arguments directly (ltmle does this for
    # SL.library, but not for nodes)
    nodes <- list(Wnodes = Wnodes, Anode = Anode, Ynode = Ynode, Vnodes = Vnodes)
    
    # fit Q and g
    message_verbose("Fitting Q", 1, verbose)
    Q_fit <- origami_SuperLearner(folds = folds, data[, Ynode], data[, c(Anode, Wnodes)], 
        family = binomial(), SL.library = SL.library$Q, cts.num = 5, .parallel = parallel, 
        method = method.NNloglik())
    Q_fit <- drop_zero_learners(Q_fit)
    
    message_verbose("Fitting g", 1, verbose)
    g_fit <- origami_SuperLearner(folds = folds, data[, Anode], data[, Wnodes], SL.library = SL.library$g, 
        family = list(family = "multinomial"), method = method.mnNNloglik())
    g_fit <- drop_zero_learners(g_fit)
    fits <- list(Q_fit = Q_fit, g_fit = g_fit)
    
    # fit rule
    message_verbose("Fitting rule", 1, verbose)
    
    # get split-specific predictions, as well as class and weight for rule learning
    message_verbose("Getting split-specific predictions", 2, verbose)
    split_preds <- cross_validate(opttx_split_preds, folds, data, nodes, fits, .combine = F, 
        .parallel = parallel)
    
    val_preds <- extract_vals(folds, split_preds)
    
    
    
    fits$rule_fit <- learn_rule(data, folds, nodes, split_preds, val_preds, parallel = F, 
        SL.library = SL.library, verbose, ...)
    
    # estimate performance
    cv_dV <- predict(fits$rule_fit, newdata = "cv-original", pred_fit = "joint")
    dV <- predict(fits$rule_fit, newdata = data[, Wnodes], pred_fit = "joint")
    message_verbose("Estimating performance", 1, verbose)
    cv_ests <- NULL
    if (perf_cv) {
        cv_ests <- estimate_performance(data, nodes, val_preds, cv_dV, perf_tmle, 
            perf_dripcw)
        cv_ests$estimator <- sprintf("CV-%s", cv_ests$estimator)
    }
    
    full_ests <- NULL
    if (perf_full) {
        predictions <- opttx_split_preds(folds[[1]], data, nodes, fits, use_full = T)
        full_ests <- estimate_performance(data, nodes, predictions, dV, perf_tmle, 
            perf_dripcw)
    }
    
    ests <- rbind(cv_ests, full_ests)
    result <- list(data = data, fits = fits, ests = ests, nodes = nodes, SL.library = SL.library, 
        folds = folds, split_preds = split_preds, val_preds = val_preds)
    result$dV <- dV
    class(result) <- "opt_tmle"
    
    return(result)
}

#' @rdname opt_tmle
#' @export
print.opt_tmle <- function(obj) {
    A <- obj$data[, obj$nodes$Anode]
    dV <- obj$dV
    n <- length(A)
    tab <- table(`observed tx` = A, `optimal tx` = dV)/n
    tab <- addmargins(tab)
    cat("Treatment Assignments\n\n")
    
    print(tab)
    cat("\n\n")
    ests <- obj$ests
    ests <- ests[order(ests$rule, ests$estimator), ]
    names(ests) <- c("Estimate", "SD", "CI Lower Bound", "CI Upper Bound", "Intervention", 
        "Estimator")
    cat("EYa Estimates\n\n")
    print(ests)
}


#' @rdname opt_tmle
#' @export
plot.opt_tmle <- function(obj) {
    ggplot(obj$ests, aes(y = rule, x = est, xmin = lower, xmax = upper)) + geom_point() + 
        geom_errorbarh() + theme_bw() + ylab("Rule") + xlab(expression(E ~ bgroup("[", 
        Y[A], "]"))) + facet_wrap(~estimator)
}

#' @rdname opt_tmle
#' @export
predict.opt_tmle <- function(obj, ...) {
    predict(obj$fits$rule_fit, ...)
} 
