library(boot)
############## somewhat general TMLE framework takes initial data, an estimation function and an update function can be used for
############## arbitrary TMLEs

# function to estimate logistic parametric submodel and get updated estimate logistic fluctuation
fluctuate <- function(tmledata, flucmod, subset = seq_len(nrow(tmledata))) {
    suppressWarnings({
        fluc <- glm(flucmod, data = tmledata[subset, ], family = "binomial")
    })
    list(eps = coef(fluc))
}

# iteratively perform tmle until convergence
gentmle <- function(initdata, estimate_fun, update_fun, max_iter = 100) {
    converge <- F
    
    # cat(sprintf('bw: %f\n',bw))
    eststep <- estimate_fun(initdata)
    initests <- eststep$ests
    order <- 1/nrow(initdata)
    for (j in seq_len(max_iter)) {
        
        updatestep <- update_fun(eststep$tmledata)
        eststep <- estimate_fun(updatestep$tmledata)
        
        ED <- sapply(eststep$Dstar, mean)
        # cat(sprintf('ED_psi=%e ED_sigma=%e psi=%f sigma2=%f\n coef_h=%f coef_Cy=%f
        # coef_Cg=%f\n',ED[1],ED[2],eststep$ests[1],sigma2=eststep$ests[2],updatestep$coefs[1],updatestep$coefs[2],updatestep$coefs[3]))
        
        
        
        if (all(abs(ED) < order)) {
            converge <- T
            break
        }
        
        
        
    }
    
    ED2 <- sapply(eststep$Dstar, function(x) mean(x^2))
    ED3 <- sapply(eststep$Dstar, function(x) mean(x^3))
    list(initdata = initdata, tmledata = eststep$tmledata, initests = initests, tmleests = eststep$ests, steps = j, Dstar = eststep$Dstar, 
        ED = ED, ED2 = ED2, ED3 = ED3)
    
    
}

############## TMLE targeting EY1

ey1_update <- function(tmledata) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < Qk & Qk < 1 & A == 1))
    eps_q <- 0
    if (length(subset) > 0) {
        # fluctuate Q
        qfluc <- fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(Qk)), subset)
        eps_q <- qfluc$eps
        tmledata$Qk <- with(tmledata, plogis(qlogis(Qk) + H1 * eps_q))
        # tmledata$Qk=qfluc$update
    }
    
    
    list(tmledata = tmledata, coefs = c(eps_q))
    
}

ey1_estimate <- function(tmledata) {
    
    psi <- mean(tmledata$Qk)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$HA <- with(tmledata, (A * H1))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - Qk) + Qk - psi)
    
    list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))
    
}

############## TMLE targeting EY1 (more canonical than above, but should be equivalent)


ey1_update2 <- function(tmledata) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < QAk & QAk < 1))
    eps_q <- 0
    if (length(subset) > 0) {
        # fluctuate Q
        qfluc <- fluctuate(tmledata, Y ~ -1 + HA + offset(qlogis(QAk)), subset)
        eps_q <- qfluc$eps
        tmledata$Q1k <- with(tmledata, plogis(qlogis(Q1k) + H1 * eps_q))
        tmledata$QAk <- with(tmledata, plogis(qlogis(QAk) + HA * eps_q))
        # tmledata$Qk=qfluc$update
    }
    
    
    list(tmledata = tmledata, coefs = c(eps_q))
    
}

ey1_estimate2 <- function(tmledata) {
    
    psi <- mean(tmledata$Q1k)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$HA <- with(tmledata, (A * H1))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - QAk) + Q1k - psi)
    
    list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))
    
}


############## TMLE targeting Sigma^2(Ey1)

sigma_update <- function(tmledata) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < Qk & Qk < 1 & A == 1))
    eps_q <- 0
    if (length(subset) > 0) {
        # fluctuate Q
        qfluc <- fluctuate(tmledata, Y ~ -1 + CyA + offset(qlogis(Qk)), subset)
        eps_q <- qfluc$eps
        tmledata$Qk <- with(tmledata, plogis(qlogis(Qk) + Cy1 * eps_q))
    }
    
    # fluctuate g tmledata$Ca=with(tmledata,Qk*(1-Qk)/(gk^2))
    gfluc <- fluctuate(tmledata, A ~ -1 + Ca + offset(qlogis(gk)))
    eps_g <- gfluc$eps
    tmledata$gk <- with(tmledata, plogis(qlogis(gk) + Ca * eps_g))
    
    list(tmledata = tmledata, coefs = c(eps_q, eps_g))
    
}

sigma_estimate <- function(tmledata) {
    
    psi <- mean(tmledata$Qk)
    
    tmledata$Cy1 <- with(tmledata, (1/gk) * ((1 - 2 * Qk)/gk + 2 * (Qk - psi)))
    tmledata$CyA <- with(tmledata, A * Cy1)
    tmledata$Ca <- with(tmledata, Qk * (1 - Qk)/(gk^2))
    
    # influence curve
    varterm <- with(tmledata, Qk * (1 - Qk)/gk + (Qk - psi)^2)
    sigma2 <- mean(varterm)
    D_Qw <- varterm - sigma2
    D_Qbar <- with(tmledata, CyA * (Y - Qk))
    D_gbar <- with(tmledata, -Ca * (A - gk))
    Dstar_sigma <- D_Qw + D_Qbar + D_gbar
    
    
    list(tmledata = tmledata, ests = c(sigma2 = sigma2), Dstar = list(Dstar_sigma = Dstar_sigma))
}

############## TMLE targeting Ey1 and Sigma^2(Ey1) simultaneously
eysigma_update <- function(tmledata) {
    # fix points where Q is already 0 or 1 - perfect prediction
    subset <- with(tmledata, which(0 < Qk & Qk < 1 & A == 1))
    eps_q <- 0
    if (length(subset) > 0) {
        # fluctuate Q
        qfluc <- fluctuate(tmledata, Y ~ -1 + HA + CyA + offset(qlogis(Qk)), subset)
        eps_q <- qfluc$eps
        tmledata$Qk <- with(tmledata, plogis(qlogis(Qk) + H1 * eps_q[1] + Cy1 * eps_q[2]))
    }
    
    # fluctuate g tmledata$Ca=with(tmledata,Qk*(1-Qk)/(gk^2))
    gfluc <- fluctuate(tmledata, A ~ -1 + Ca + offset(qlogis(gk)))
    eps_g <- gfluc$eps
    tmledata$gk <- with(tmledata, plogis(qlogis(gk) + Ca * eps_g))
    
    list(tmledata = tmledata, coefs = c(eps_q, eps_g))
    
}

eysigma_estimate <- function(tmledata) {
    
    psi <- mean(tmledata$Qk)
    
    tmledata$H1 <- with(tmledata, (1/gk))
    tmledata$HA <- with(tmledata, (A * H1))
    
    tmledata$Cy1 <- with(tmledata, (1/gk) * ((1 - 2 * Qk)/gk + 2 * (Qk - psi)))
    tmledata$CyA <- with(tmledata, A * Cy1)
    
    tmledata$Ca <- with(tmledata, Qk * (1 - Qk)/(gk^2))
    
    # influence curves
    Dstar_psi <- with(tmledata, HA * (Y - Qk) + Qk - psi)
    
    varterm <- with(tmledata, Qk * (1 - Qk)/gk + (Qk - psi)^2)
    sigma2 <- mean(varterm)
    D_Qw <- varterm - sigma2
    D_Qbar <- with(tmledata, CyA * (Y - Qk))
    D_gbar <- with(tmledata, -Ca * (A - gk))
    Dstar_sigma <- D_Qw + D_Qbar + D_gbar
    
    
    
    list(tmledata = tmledata, ests = c(psi = psi, sigma2 = sigma2), Dstar = list(Dstar_psi = Dstar_psi, Dstar_sigma = Dstar_sigma))
}


# example code to construct initial tmle data A=data$a Y=data$y W=data$w g1=estg(W,A)#g(data$wstar)
# Q1=ksmooth(knnfit,bw) order=1/length(Y) Qk=Q1 gk=g1 psi=mean(Qk) var1=mean(Qk*(1-Qk)/gk+(Qk-psi)^2)
# tmledata=data.frame(Qk,gk,A,Y)
TRUE 
