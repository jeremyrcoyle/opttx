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
        folds <- make_folds(V = 10)
    }
    
    # Y=data[,nodes$Ynode] X=data[,c(nodes$Anode,nodes$Wnodes)] osl_args<-
    # list(family = binomial(), SL.library = SL.library$Q, cts.num = 5, .parallel =
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
        nodes$Wnodes)], family = binomial(), SL.library = SL.library$Q, cts.num = 5, 
        .parallel = parallel, method = method.NNloglik(), control = list(trimLogit = 1e-05))
    Q_fit <- split_from_args(Q_fit_args)
    # Q_fit_full <- split_to_full(Q_fit) Q_fit_nested <-
    # split_to_nested(Q_fit,Q_fit_args)
    #Q_fit <- drop_zero_learners(Q_fit)
    
    message_verbose("Fitting g", 1, verbose)
    g_fit_args <- list(folds = folds, Y = data[, nodes$Anode], X = data[, nodes$Wnodes], 
        SL.library = SL.library$g, family = list(family = "multinomial"), method = method.mnNNloglik())
    g_fit <- split_from_args(g_fit_args)
    #g_fit <- drop_zero_learners(g_fit)
    fits <- list(Q_fit = Q_fit, g_fit = g_fit)
    
    # fit rule
    message_verbose("Fitting rule", 1, verbose)
    
    # get split-specific predictions, as well as class and weight for rule learning
    message_verbose("Getting split-specific predictions", 2, verbose)
    split_preds <- cross_validate(opttx_split_preds, folds, data, nodes, fits, .combine = F, 
        .parallel = parallel)
    
    full_preds <- opttx_split_preds(folds[[1]], data, nodes, fits, use_full = T)
    val_preds <- extract_vals(folds, split_preds)
    
    split_folds=make_folds(strata_ids=val_preds$v)
    
    split_preds2 <- cross_validate(opttx_split_preds, split_folds, data, nodes, fits, .combine = F, 
        .parallel = parallel)
    val_preds2 <- extract_vals(split_folds, split_preds2)
    
    table(val_preds$v,val_preds2$v)
    
    # fits$rule_fit <- learn_rule(data, folds, nodes, split_preds, full_preds,
    # val_preds, parallel = F, SL.library = SL.library, verbose) #, ...)
    nuisance_preds <- val_preds
    nuisance_preds$A <- data[, Anode]
    nuisance_preds$Y <- data[, Ynode]
    
    method <- method.EYd(nuisance_preds)
    QaV_fit <- origami_SuperLearner(folds = folds, data[, nodes$Ynode], data[, nodes$Vnodes, 
        drop = F], split_preds = split_preds, full_preds = full_preds, SL.library = SL.library$QaV, 
        family = gaussian(), cvfun = QaV_cv_SL, .parallel = parallel, method = method)

    QaV_fit2 <- origami_SuperLearner(folds = folds, data[, nodes$Ynode], data[, nodes$Vnodes, 
        drop = F], split_preds = split_preds, full_preds = full_preds, SL.library = SL.library$QaV, 
        family = gaussian(), cvfun = QaV_cv_SL, .parallel = parallel, method = method.mvSL(method.NNLS()))

    # estimate performance cv_dV <- predict(fits$rule_fit, newdata = 'cv-original',
    # pred_fit = 'joint')
    cv_dV <- max.col(predict(QaV_fit, newdata = "cv-original")$pred)
    
    split_folds <- make_folds(data, V = 10)
    
    newQaV <- cross_validate(refit_split, split_folds, QaV_fit)
    QaV_split <- index_dim(newQaV$preds, order(newQaV$index))
    split_fold=split_folds[[1]]
    test=function(split_fold,folds){
        newQ <- refit_split(split_fold,fits$Q_fit, control = fits$Q_fit$fullFit$control)$fit[[1]]
        newg <- refit_split(split_fold,fits$g_fit, control = fits$g_fit$fullFit$control)$fit[[1]]
        new_fits=list(Q_fit=newQ,g_fit=newg)
        valid_data=validation(data)
        fold=folds[[1]]
        split_preds2 <- cross_validate(opttx_split_preds, folds, data, nodes, new_fits, .combine = F, 
        .parallel = parallel)
        
        QaV_fit2 <- origami_SuperLearner(folds = folds, data[, nodes$Ynode], data[, nodes$Vnodes, 
        drop = F], split_preds = split_preds2, full_preds = full_preds, SL.library = SL.library$QaV, 
        family = gaussian(), cvfun = QaV_cv_SL, .parallel = parallel, method = method)
        
        #refit_weights
        
        train_preds=lapply(split_preds2,function(x){lapply(x,training,fold=split_fold)})
        str(train_preds)
        method2=method.EYd(train_preds)
        QaV_fit$fullFit$method=method2
        newQaV <- refit_split(split_fold, QaV_fit, control = fits$Q_fit$fullFit$control)
        str(split_preds2)
        validation(fold=split_fold)%in%validation(fold=folds[[2]])
    }
    newg <- cross_validate(refit_split, split_folds, fits$g_fit, control = fits$g_fit$fullFit$control)
    g_split <- index_dim(newg$preds, order(newg$index))
    
    cv_dV2 <- max.col(QaV_split)
    
    
    
    test_dV <- max.col(predict(QaV_fit, newdata = testdata[, Wnodes])$pred)
    mean(Qbar0(test_dV, testdata[, Wnodes]))
    
    test_dV2 <- max.col(predict(QaV_fit2, newdata = testdata[, Wnodes])$pred)
    mean(Qbar0(test_dV2, testdata[, Wnodes]))
    
    table(eyd=test_dV==testdata$d0, nnls=test_dV2==testdata$d0)
    
    mean(Qbar0(testdata$d0, testdata[, Wnodes]))
    rule_tmle(data$A, data$Y, val_preds$pA, val_preds$QaW, cv_dV)
    rule_tmle(data$A, data$Y, val_preds$pA, val_preds$QaW, cv_dV2)
    rule_tmle(data$A, data$Y, g_split, val_preds$QaW, cv_dV2)
    dV <- predict(fits$rule_fit, newdata = data[, nodes$Wnodes], pred_fit = "joint")
    message_verbose("Estimating performance", 1, verbose)
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
