
# split specific predictions for Q and g and derived quantities
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
    
    pA <- predict(splitg_fit, data[, nodes$Wnodes])$pred
    
    # split specific blip, class, and weights
    A <- data[, nodes$Anode]
    Y <- data[, nodes$Ynode]
    A_ind <- factor_to_indicators(A)
    Y_mat <- replicate(length(A_vals), Y)
    DR <- (A_ind/pA) * (Y_mat - QaW) + QaW
    
    Z <- max.col(DR)
    
    list(QaW = QaW, pA = pA, DR = DR, Z = Z)  #, Y=Y, A=A)
}

# ignores Y and instead uses Z generated from a split-specific training set
# generated using cv_split_preds X are the nodes to base the rule on
#' @export
QaV_cv_SL <- function(fold, Y, X, SL.library, family, obsWeights, id, use_full = F, 
    split_preds, A_index, blip_type = "DR", ...) {
    v <- fold_index()
    train_idx <- training()
    valid_idx <- validation()
    
    DR <- split_preds[["DR"]][[v]]
    DR_a <- DR[, A_index]
    if (blip_type == "blip1") {
        to_predict <- DR_a - DR[, 1]
    } else if (blip_type == "blip2") {
        to_predict <- DR_a - rowMeans(DR)
    } else if (blip_type == "blip3") {
        pA <- split_preds[["pA"]][[v]]
        to_predict <- DR_a - rowSums(DR * pA)
    } else {
        to_predict <- DR_a
    }
    
    cv_SL(fold, to_predict, X, SL.library, family, obsWeights, id)
}

#' @export
class_cv_SL <- function(fold, Y, X, SL.library, family, obsWeights, id, use_full = F, 
    split_preds, ...) {
    v <- fold_index()
    train_idx <- training()
    valid_idx <- validation()
    
    Z <- split_preds$Z[[v]]
    
    cv_SL(fold, Z, X, SL.library, family, obsWeights, id)
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
