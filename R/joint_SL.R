

create_predmat <- function(QaV_fit, class_fit, newdata = "cv-original") {
    QaV_Z <- NULL
    QaV_coef <- c()
    if (!is.null(QaV_fit)) {
        QaV_Z <- predict(QaV_fit, newdata)$library_pred
        
        QaV_coef <- QaV_fit$coef
    }
    
    class_Z <- NULL
    class_coef <- c()
    if (!is.null(class_fit)) {
        class_Z <- predict(class_fit, newdata)$library_pred
        # transform class_Z
        class_Z <- trimLogit(class_Z)
        class_coef <- class_fit$coef
    }
    
    combined_Z <- abind(QaV_Z, class_Z, rev.along = 1)
    
    QaV_num_alg <- length(QaV_coef)
    class_num_alg <- length(class_coef)
    combined_num_alg <- QaV_num_alg + class_num_alg
    
    combined_coef <- c((QaV_num_alg/combined_num_alg) * QaV_coef, (class_num_alg/combined_num_alg) * 
        class_coef)
    
    list(Z = combined_Z, init_coef = combined_coef)
    
}

# combine QaV and class_fits
joint_sl <- function(QaV_fit, class_fit, predictions, data, nodes, risk_generator = create_tmle_risk) {
    
    # evaluate combined coefficient, plus a small number of random starting points to
    # find a good neighborhood
    jsl_obj <- create_predmat(QaV_fit, class_fit, "cv-original")
    num_alg <- length(jsl_obj$init_coef)
    risk_fun <- risk_generator(jsl_obj$Z, data, nodes, predictions)
    simplex.grid <- simplex.sample(num_alg, 30)$samples
    starts <- rbind(simplex.grid, jsl_obj$init_coef)
    start_risk <- apply(starts, 1, risk_fun)
    optim_init <- starts[which.min(start_risk), ]
    jsl_obj$grid_coef <- optim_init
    
    # optimize starting in that neighborhood
    n <- dim(jsl_obj$Z)[1]
    fit <- nloptr(x0 = optim_init, eval_f = risk_fun, lb = rep(0, num_alg), opts = list(algorithm = "NLOPT_LN_SBPLX", 
        ftol_rel = 1/n, maxeval = n))
    
    jsl_obj$coef <- normalize(fit$solution)
    class(jsl_obj) <- "joint_sl"
    return(jsl_obj)
}

predict.joint_sl <- function(jsl_obj, QaV_fits, class_fit, newdata = "cv-original") {
    Z <- create_predmat(QaV_fits, class_fit, newdata)$Z
    
    mn_pred(jsl_obj$coef, Z)
} 
