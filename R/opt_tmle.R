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
opt_tmle <- function(data, Wnodes = grep("^W", names(data), value = TRUE), Anode = "A", 
    Ynode = "Y", Vnodes = Wnodes, stratifyAY = TRUE, SL.library = opt_tmle.SL.library, 
    verbose = 3, parallel = FALSE, perf_tmle = TRUE, perf_dripcw = FALSE, perf_cv = TRUE, 
    perf_full = FALSE, maximize = TRUE, ...) {
    
    # ensure A is a factor
    data[, Anode] <- as.factor(data[, Anode])
    
    if (stratifyAY) {
        AYstrata <- sprintf("%s %s", data[, Anode], data[, Ynode])
        folds <- make_folds(strata_ids = AYstrata, V = 10)
    } else {
        folds <- make_folds(data, V = 10)
    }
    
    # Y=data[,nodes$Ynode] X=data[,c(nodes$Anode,nodes$Wnodes)] osl_args<-
    # list(family = binomial(), SL.library = SL.library$Q, cts.num = 10, .parallel =
    # parallel, method = method.NNloglik(), control = list(trimLogit = 1e-05))
    
    # split=split_SL(folds,Y,X,osl_args)
    # full=split_SL(folds,Y,X,osl_args,approach='full')
    # nested=split_SL(folds,Y,X,osl_args,approach='full')
    
    data$Ystar <- data[, Ynode]
    if (!maximize) {
        minY <- min(data$Ystar)
        maxY <- max(data$Ystar)
        data$Ystar <- (minY + maxY) - data$Ystar
    }
    
    # possibly we should make these lists the arguments directly (ltmle does this for
    # SL.library, but not for nodes)
    nodes <- list(Wnodes = Wnodes, Anode = Anode, Ynode = "Ystar", Vnodes = Vnodes)
    
    # fit Q and g
    message_verbose("Fitting Q", 1, verbose)
    # todo: add support for continuous Y
    Q_fit_args <- list(folds = folds, Y = data[, nodes$Ynode], X = data[, c(nodes$Anode, 
        nodes$Wnodes)], family = binomial(), SL.library = SL.library$Q, cts.num = 10, 
        .parallel = parallel, method = method.NNloglik(), control = list(trimLogit = 1e-05))
    Q_fit <- split_from_args(Q_fit_args)
    
    # Q_fit_full <- split_to_full(Q_fit) Q_fit_nested <-
    # split_to_nested(Q_fit,Q_fit_args) Q_fit <- drop_zero_learners(Q_fit)
    
    message_verbose("Fitting g", 1, verbose)
    g_fit_args <- list(folds = folds, Y = data[, nodes$Anode], X = data[, nodes$Wnodes], 
        SL.library = SL.library$g, family = list(family = "multinomial"), method = method.mnNNloglik())
    g_fit <- split_from_args(g_fit_args)
    # g_fit <- drop_zero_learners(g_fit)
    fits <- list(Q_fit = Q_fit, g_fit = g_fit)
    
    # get split-specific predictions, as well as class and weight for rule learning
    message_verbose("Getting split-specific predictions", 2, verbose)
    split_preds <- cross_validate(opttx_split_preds, folds, data, nodes, fits, .combine = F, 
        .parallel = parallel)
    
    full_preds <- opttx_split_preds(folds[[1]], data, nodes, fits, use_full = T)
    val_preds <- extract_vals(folds, split_preds)
    EYdmethod <- method.EYd(val_preds)
    mvSLmethod <- method.mvSL(method.NNLS())
    # debug(mvSLmethod$computeCoef) fit rule
    message_verbose("Fitting rule", 1, verbose)
    # full_preds$DR=(plogis(trimLogit(full_preds$DR)))
    rule_args <- list(folds = folds, Y = data[, nodes$Ynode], X = data[, nodes$Vnodes, 
        drop = F], cts.num = 10, split_preds = split_preds, full_preds = full_preds, 
        SL.library = SL.library$QaV, family = gaussian(), cvfun = QaV_cv_SL, .parallel = parallel, 
        blip_type = "DR", use_full = F, method = mvSLmethod)  # EYdmethod)   # mvSLmethod)  # 
    
    fits$rule_fit <- split_from_args(rule_args)
    
    # estimate performance
    message_verbose("Estimating performance", 1, verbose)
    A_vals <- vals_from_factor(data[, nodes$Anode])
    cv_dV <- dV_from_preds(predict(fits$rule_fit, newdata = "cv-original")$pred, 
        A_vals)
    dV <- dV_from_preds(predict(fits$rule_fit, newdata = data[, nodes$Vnodes])$pred, 
        A_vals)
    
    
    cv_ests <- NULL
    nodes$Ynode <- Ynode
    if (perf_cv) {
        cv_ests <- estimate_performance(data, nodes, val_preds, cv_dV, perf_tmle, 
            perf_dripcw)
        cv_ests$estimator <- sprintf("CV-%s", cv_ests$estimator)
    }
    
    full_ests <- NULL
    if (perf_full) {
        full_ests <- estimate_performance(data, nodes, full_preds, dV, perf_tmle, 
            perf_dripcw)
    }
    
    ests <- rbind(cv_ests, full_ests)
    result <- list(data = data, fits = fits, ests = ests, nodes = nodes, SL.library = SL.library, 
        folds = folds, split_preds = split_preds, val_preds = val_preds, full_preds = full_preds, 
        rule_args = rule_args, dV = dV, cv_dV = cv_dV)
    
    class(result) <- "opt_tmle"
    plot(result)
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
    
    print(tab, width = 200)
    cat("\n\n")
    ests <- obj$ests
    ests <- ests[order(ests$estimator, ests$rule), ]
    names(ests) <- c("Estimate", "SE", "CI Lower Bound", "CI Upper Bound", "Intervention", 
        "Estimator")
    cat("EYa Estimates\n\n")
    print(ests, width = 200)
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
    pred_QaV <- predict(obj$fits$rule_fit, ...)$pred
    
    dV <- dV_from_preds(pred_QaV)
    
    list(pred_QaV = pred_QaV, dV = dV)
} 
