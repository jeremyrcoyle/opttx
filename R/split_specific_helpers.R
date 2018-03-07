
# split specific predictions for Q and g and derived quantities
#' export
opttx_split_preds <- function(fold, data, nodes, fits, use_full = F, ...) {
    v <- fold_index()
    train_idx <- training()
    valid_idx <- validation()
    
    if ((all.equal(train_idx, valid_idx) == T) || use_full) {
        # we're in final resubstitution call, so let's use full Q and g fits
        splitQ_fit <- fits$Q_fit$fullFit
        splitg_fit <- fits$g_fit$fullFit
    } else {
        # split specific Super Learners
        splitQ_fit <- fits$Q_fit$foldFits[[v]]
        splitg_fit <- fits$g_fit$foldFits[[v]]
    }
    
    # split specific estimates
    A_vals <- vals_from_factor(data[, nodes$Anode])
    QaW <- pred_all_Q(splitQ_fit, data[, c(nodes$Anode, nodes$Wnodes)], A_vals, nodes$Anode)
    newdata <- data[, c(nodes$Anode, nodes$Wnodes)]
    pA <- predict(splitg_fit, data[, nodes$Wnodes, drop = F])$pred
    pA[pA < 0.05] <- 0.05
    if (ncol(pA) == 1) {
        pA <- cbind(1 - pA, pA)
    }
    # split specific blip, class, and weights
    A <- data[, nodes$Anode]
    Y <- data[, nodes$Ynode]
    A_ind <- factor_to_indicators(A)
    Y_mat <- replicate(length(A_vals), Y)
    DR <- (A_ind/pA) * (Y_mat - QaW) + QaW
    
    Z <- max.col(DR)
    
    list(QaW = QaW, pA = pA, DR = DR, Z = Z, Y = Y, A = A, v = rep(v, length(Y)))
}

opttx_refit_split_preds <- function(fold, data, nodes, fits, use_full = F, ...) {
    v <- fold_index()
    train_idx <- training()
    valid_idx <- validation()
    
    if ((all.equal(train_idx, valid_idx) == T) || use_full) {
        # we're in final resubstitution call, so let's use full Q and g fits
        splitrule_fit <- fits$rule_fit$fullFit
    } else {
        # split specific Super Learners
        splitrule_fit <- fits$rule_fit$foldFits[[v]]
    }
    
    # split specific estimates
    DR <- predict(splitrule_fit, newdata = (data[, nodes$Vnodes]))$pred
    Z <- max.col(DR)
    
    list(DR = DR, Z = Z, v = rep(v, length(Z)))
}


# extract the validation sets from split_preds for use in CV-TMLE
extract_val <- function(fold, split_preds) {
    v <- fold_index()
    valid_idx <- validation()
    
    val_preds <- sapply(split_preds, function(split_pred) {
        index_dim(split_pred[[v]], valid_idx)
    })
    val_preds$index <- valid_idx
    
    return(val_preds)
}

# extract split-specific validation predictions
extract_vals <- function(folds, split_preds, parallel = F) {
    
    val_preds <- cross_validate(extract_val, folds, split_preds, .parallel = parallel)
    
    data_order <- order(val_preds$index)
    
    lapply(val_preds, index_dim, data_order)
    
}
# todo: reintroduce binary versions of these based on blips
QaV_cv_SL <- function(fold, Y, X, SL.library, family, obsWeights, id, split_preds, full_preds, blip_type = "DR", 
    use_full = F, ...) {
    v <- fold_index()
    train_idx <- training()
    valid_idx <- validation()
    
    if ((all.equal(train_idx, valid_idx) == T) || use_full) {
        preds <- full_preds
    } else {
        preds <- lapply(split_preds, "[[", v)
    }
    
    DR <- preds[["DR"]]
    if (blip_type == "blip1") {
        to_predict <- DR[, -1] - DR[, 1]
    } else if (blip_type == "blip2") {
        to_predict <- DR - rowMeans(DR)
    } else if (blip_type == "blip3") {
        pA <- preds[["pA"]]
        to_predict <- DR - rowSums(DR * pA)
    } else {
        to_predict <- DR
    }
    
    cv_SL(fold, to_predict, X, SL.library, family, obsWeights, id, ...)
}
#' @export
class_cv_SL <- function(fold, Y, X, SL.library, family, obsWeights, id, split_preds, full_preds, use_full = F, 
    ...) {
    v <- fold_index()
    train_idx <- training()
    valid_idx <- validation()
    
    
    if ((all.equal(train_idx, valid_idx) == T) || use_full) {
        preds <- full_preds
    } else {
        preds <- lapply(split_preds, "[[", v)
    }
    
    Z <- preds$Z
    
    cv_SL(fold, Z, X, SL.library, family, obsWeights, id)
} 
