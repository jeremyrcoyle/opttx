library(opttx)
library(devtools)
load_all("~/opttx")
setwd("~/opttx/inst/simulations")
Qbar0 <- function(A, W) {
    
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    W4 <- W[, 4]
    Qbar <- plogis((A == 1) * (1 + 5 * W1^2 - 4.45) + (A == 2) * (-1 + 4 * W2) + 
        (1 + 5 * W3) + 5 * W3 * sin(W4))
    return(Qbar)
}

g0 <- function(W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    W4 <- W[, 4]
    
    # rep(0.5, nrow(W))
    scale_factor <- 0.8
    A1 <- plogis(scale_factor * W1)
    A2 <- plogis(scale_factor * W2)
    A3 <- plogis(scale_factor * W3)
    A <- cbind(A1, A2, A3)
    
    # make sure A sums to 1
    A <- normalize_rows(A)
}

gen_data <- function(n = 1000, p = 4) {
    W <- matrix(rnorm(n * p), nrow = n)
    colnames(W) <- paste("W", seq_len(p), sep = "")
    pA <- g0(W)
    A <- factor(apply(pA, 1, function(pAi) which(rmultinom(1, 1, pAi) == 1)))
    A_vals <- vals_from_factor(A)
    
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(A, W))
    Q0aW <- sapply(A_vals, Qbar0, W)
    d0 <- max.col(Q0aW)
    Yd0 <- as.numeric(u < Qbar0(d0, W))
    data.frame(W, A, Y, Q0aW, d0, Yd0)
}





SL.library <- list(Q = c("SL.glm", "SL.glmem", "SL.glmnetprob", "SL.step.forward", 
    "SL.gam", "SL.rpart", "SL.rpartPrune", "SL.mean"), g = c("mnSL.randomForest", 
    "mnSL.glmnet", "mnSL.multinom", "mnSL.mean"), QaV = c("SL.glm", "SL.glmnetprob", 
    "SL.step.forward", "SL.gam", "SL.rpart", "SL.rpartPrune", "SL.mean"), class = c("mnSL.randomForest", 
    "mnSL.glmnet", "mnSL.multinom", "mnSL.mean"))

SL.library$QaV <- sl_to_mv_library(SL.library$QaV)
parallel <- FALSE
verbose <- 4


nestcvtmle <- function(fold, data, nodes) {
    traindata <- training(data)
    valdata <- validation(data)
    
    AYstrata <- sprintf("%s %s", traindata$A, traindata$Y)
    nestfolds <- make_folds(strata_ids = AYstrata, V = 10)
    
    Q_fit_args <- list(folds = nestfolds, Y = traindata[, nodes$Ynode], X = traindata[, 
        c(nodes$Anode, nodes$Wnodes)], family = binomial(), SL.library = SL.library$Q, 
        cts.num = 5, .parallel = parallel, method = method.NNloglik(), control = list(trimLogit = 1e-05))
    Q_fit <- split_from_args(Q_fit_args)
    Q_fit <- drop_zero_learners(Q_fit)
    
    g_fit_args <- list(folds = nestfolds, Y = traindata[, nodes$Anode], X = traindata[, 
        nodes$Wnodes], SL.library = SL.library$g, family = list(family = "multinomial"), 
        method = method.mnNNloglik())
    g_fit <- split_from_args(g_fit_args)
    g_fit <- drop_zero_learners(g_fit)
    fits <- list(Q_fit = Q_fit, g_fit = g_fit)
    
    split_preds <- cross_validate(opttx_split_preds, nestfolds, traindata, nodes, 
        fits, .combine = F, .parallel = parallel)
    
    full_preds <- opttx_split_preds(nestfolds[[1]], traindata, nodes, fits, use_full = T)
    val_preds <- extract_vals(nestfolds, split_preds)
    
    fits$split_rule_fit <- learn_rule(traindata, nestfolds, nodes, split_preds, full_preds, 
        val_preds, parallel = F, SL.library = SL.library, verbose)
    
    fits$full_rule_fit <- learn_rule(traindata, nestfolds, nodes, split_preds, full_preds, 
        val_preds, parallel = F, SL.library = SL.library, verbose, use_full = T)
    
    outer_val_preds <- opttx_split_preds(nestfolds[[1]], valdata, nodes, fits, use_full = T)
    
    outer_val_split_dV <- predict(fits$split_rule_fit, newdata = valdata[, nodes$Vnodes])
    outer_val_full_dV <- predict(fits$full_rule_fit, newdata = valdata[, nodes$Vnodes])
    c(list(A = valdata$A, Y = valdata$Y, split_dV = outer_val_split_dV, full_dV = outer_val_full_dV), 
        outer_val_preds)
}

n = 1000
p = 5
testdata <- gen_data(1e+06, p)

Wnodes <- grep("^W", names(testdata), value = TRUE)
Anode <- "A"
Ynode <- "Y"
Vnodes <- Wnodes
nodes <- list(Wnodes = Wnodes, Anode = Anode, Ynode = Ynode, Vnodes = Vnodes)

test_EYd0 <- mean(testdata$Yd0)
QaV0 <- sapply(unique(testdata$A), Qbar0, testdata[, nodes$Vnode])

sim <- function(iteration, n = 1000, p = 5) {
    
    
    data <- gen_data(n, p)
    
    data[, Anode] <- as.factor(data[, Anode])
    
    AYstrata <- sprintf("%s %s", data$A, data$Y)
    folds <- make_folds(strata_ids = AYstrata, V = 10)
    
    cat(sprintf("%d: Q and g fits\n", iteration))
    # fit Q and g on full data
    qgtime <- system.time({
        
        # fit Q and g
        message_verbose("Fitting Q", 1, verbose)
        # todo: add support for continuous Y
        Q_fit_args <- list(folds = folds, Y = data[, nodes$Ynode], X = data[, c(nodes$Anode, 
            nodes$Wnodes)], family = binomial(), SL.library = SL.library$Q, cts.num = 5, 
            .parallel = parallel, method = method.NNloglik(), control = list(trimLogit = 1e-05))
        Q_fit <- split_from_args(Q_fit_args)
        Q_fit <- drop_zero_learners(Q_fit)
        
        message_verbose("Fitting g", 1, verbose)
        g_fit_args <- list(folds = folds, Y = data[, nodes$Anode], X = data[, nodes$Wnodes], 
            SL.library = SL.library$g, family = list(family = "multinomial"), method = method.mnNNloglik())
        g_fit <- split_from_args(g_fit_args)
        g_fit <- drop_zero_learners(g_fit)
        fits <- list(Q_fit = Q_fit, g_fit = g_fit)
    })
    
    # learn rule with various methods
    cat(sprintf("%d: split predictions\n", iteration))
    splitpredtime <- system.time({
        split_preds <- cross_validate(opttx_split_preds, folds, data, nodes, fits, 
            .combine = F, .parallel = parallel)
        
        full_preds <- opttx_split_preds(folds[[1]], data, nodes, fits, use_full = T)
        val_preds <- extract_vals(folds, split_preds)
    })
    
    # split sequential
    cat(sprintf("%d: split rule\n", iteration))
    
    splitruletime <- qgtime + splitpredtime + system.time({
        fits$split_rule_fit <- learn_rule(data, folds, nodes, split_preds, full_preds, 
            val_preds, parallel = F, SL.library = SL.library, verbose)
    })
    
    # full sequential
    cat(sprintf("%d: full rule\n", iteration))
    fullruletime <- qgtime + splitpredtime + system.time({
        fits$full_rule_fit <- learn_rule(data, folds, nodes, split_preds, full_preds, 
            val_preds, parallel = F, SL.library = SL.library, verbose, use_full = T)
    })
    
    # nested
    cat(sprintf("%d: nested rule\n", iteration))
    nestruletime <- qgtime + system.time({
        
        nested_Q <- split_to_nested(Q_fit, Q_fit_args)
        nested_g <- split_to_nested(g_fit, g_fit_args)
        nested_fits <- list(Q_fit = nested_Q, g_fit = nested_g)
        
        nested_split_preds <- cross_validate(opttx_split_preds, folds, data, nodes, 
            nested_fits, .combine = F, .parallel = parallel)
        
        nested_full_preds <- opttx_split_preds(folds[[1]], data, nodes, nested_fits, 
            use_full = T)
        nested_val_preds <- extract_vals(folds, nested_split_preds)
        
        fits$nested_rule_fit <- learn_rule(data, folds, nodes, nested_split_preds, nested_full_preds, 
            nested_val_preds, parallel = F, SL.library = SL.library, verbose)
    })
    
    # test set risk
    cat(sprintf("%d: risks\n", iteration))
    
    test_split_dV <- predict(fits$split_rule_fit, newdata = testdata[, Vnodes])
    test_split_EYdn <- mean(Qbar0(test_split_dV, testdata[, Wnodes]))
    
    test_full_dV <- predict(fits$full_rule_fit, newdata = testdata[, Vnodes])
    test_full_EYdn <- mean(Qbar0(test_full_dV, testdata[, Wnodes]))
    
    test_nested_dV <- predict(fits$nested_rule_fit, newdata = testdata[, Vnodes])
    test_nested_EYdn <- mean(Qbar0(test_nested_dV, testdata[, Wnodes]))
        
    test_split_QaV <- predict(fits$split_rule_fit, newdata = testdata[, Vnodes], 
        pred_fit = "QaV", return_assignment = F)
    test_split_QaV_mse <- mean((test_split_QaV - QaV0)^2)
    
    test_full_QaV <- predict(fits$full_rule_fit, newdata = testdata[, Vnodes], pred_fit = "QaV", 
        return_assignment = F)
    test_full_QaV_mse <- mean((test_full_QaV - QaV0)^2)
    
    test_nested_QaV <- predict(fits$nested_rule_fit, newdata = testdata[, Vnodes], 
        pred_fit = "QaV", return_assignment = F)
    test_nested_QaV_mse <- mean((test_nested_QaV - QaV0)^2)
    
    test_perf <- data.frame(test_split_EYdn, test_full_EYdn, test_nested_EYdn, test_EYd0, test_split_QaV_mse,test_full_QaV_mse,test_nested_QaV_mse)
    # estimate EYdn with TMLE
    cat(sprintf("%d: tmle\n", iteration))
    
    tmletime <- qgtime + splitpredtime + system.time({
        split_dV <- predict(fits$split_rule_fit, newdata = data[, Vnodes])
        full_dV <- predict(fits$full_rule_fit, newdata = data[, Vnodes])
        nested_dV <- predict(fits$nested_rule_fit, newdata = data[, Vnodes])
        full_perf <- ldply(list(split_dV, full_dV, nested_dV), function(dV) {
            estimate_performance(data, nodes, full_preds, dV, perf_tmle = T, perf_dripcw = T)
        })
        full_perf$rule_fit <- rep(c("split", "full", "nested"), each = 10)
    })
    
    # estimate EYdn with CV-TMLE avoiding nesting (SplitSequential)
    cat(sprintf("%d: cv-tmle\n", iteration))
    
    cvtmletime <- qgtime + splitpredtime + system.time({
        cv_split_dV <- predict(fits$split_rule_fit, newdata = "cv-original")
        cv_full_dV <- predict(fits$full_rule_fit, newdata = "cv-original")
        cv_nested_dV <- predict(fits$nested_rule_fit, newdata = "cv-original")
        cv_perf <- ldply(list(cv_split_dV, cv_full_dV, cv_nested_dV), function(dV) {
            estimate_performance(data, nodes, val_preds, dV, perf_tmle = T, perf_dripcw = T)
        })
        cv_perf$rule_fit <- rep(c("split", "full", "nested"), each = 10)
        cv_perf$estimator <- sprintf("CV-%s", cv_perf$estimator)
    })
    
    # estimate EYdn with CV-TMLE with nesting
    cat(sprintf("%d: nested cv-tmle\n", iteration))
    
    nestcvtmletime <- system.time({
        cvres <- cross_validate(nestcvtmle, folds, data = data, nodes = nodes)
        nested_cv_perf <- ldply(list(cvres$full_dV, cvres$split_dV), function(dV) {
            estimate_performance(data = as.data.frame(cvres), nodes, cvres, dV, perf_tmle = T, 
                perf_dripcw = T)
        })
        nested_cv_perf$rule_fit <- rep(c("split", "full"), each = 10)
        nested_cv_perf$estimator <- sprintf("Nested CV-%s", nested_cv_perf$estimator)
    })
    
    perf_ests <- rbind.fill(full_perf, cv_perf, nested_cv_perf)
    times <- list(qgtime = qgtime, splitpred = splitpredtime, splitrule = splitruletime, 
        fullrule = fullruletime, nestrule = nestruletime, tmle = tmletime, cvtmle = cvtmletime, 
        nestcvtmle = nestcvtmletime)
    # gather results
    results <- list(iteration = iteration, times = times, perf_ests = perf_ests, 
        test_perf = test_perf)
    
    save(results, file = sprintf("simresults/nestingsim_%d.rdata", iteration))
    
    results
} 


library(doMC)
registerDoMC(16)
#iter=1
allresults <- foreach(iter = 1:1000, .options.multicore = list(preschedule = F),.errorhandling="remove") %dopar% {
    results <- sim(iter, n, p)
}

save(allresults, file = "simresults/nesting_sim_all.rdata")

