
Qbar0 <- function(A, W) {
    
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    W4 <- W[, 4]
    Qbar <- (1/2) * (plogis(-5 * (A == 2) * (W1 + 0.5) + 5 * (A == 3) * (W1 - 0.5)) + plogis(W2 * W3))
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
    g0W <- g0(W)
    A <- factor(apply(g0W, 1, function(pAi) which(rmultinom(1, 1, pAi) == 1)))
    A_vals <- vals_from_factor(A)
    
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(A, W))
    Q0aW <- sapply(A_vals, Qbar0, W)
    d0 <- max.col(Q0aW)
    Yd0 <- as.numeric(u < Qbar0(d0, W))
    df <- data.frame(W, A, Y, d0, Yd0)
    
    df$g0W <- g0(W)
    df$Q0aW <- Q0aW
    
    return(df)
}






SL.library <- list(Q = c("SL.glm", "SL.glmem", "SL.glmnet", "SL.glmnetem", "SL.polymars", "SL.step.forward", 
    "SL.gam", "SL.mean"), g = c("mnSL.glmnet", "mnSL.multinom", "mnSL.mean", "mnSL.polymars"), QaV = c("SL.polymars", 
    "SL.glm", "SL.glmnet", "SL.step.forward", "SL.gam", "SL.mean"))

SL.library$QaV <- sl_to_mv_library(SL.library$QaV)
# SL.library$Q <- sl_to_strat_library(SL.library$Q, 'A')

testdata <- gen_data(1e+05, 5)
data <- gen_data(1000, 5)

Wnodes <- grep("^W", names(data), value = TRUE)
Anode <- "A"
Ynode <- "Y"
Vnodes <- Wnodes
nodes <- list(Wnodes = Wnodes, Anode = Anode, Ynode = "Ystar", Vnodes = Vnodes)

system.time({
    result <- opt_tmle(data, SL.library = SL.library)
})

vimresult <- backward_vim(result, testdata, Qbar0)

ggplot(vimresult$vimdf, aes(y = Vnode, x = est, xmin = lower, xmax = upper)) + geom_point() + geom_point(aes(x = test), 
    color = "red") + geom_errorbarh() + facet_wrap(~metric, scales = "free") + theme_bw()

print(result)
plot(result)
Wnodes <- result$nodes$Wnodes
QaV_dV <- predict(result, newdata = testdata[, Wnodes], pred_fit = "QaV")$dV
QaV_perf <- mean(Qbar0(QaV_dV, testdata[, Wnodes]))
class_dV <- predict(result, newdata = testdata[, Wnodes], pred_fit = "class")$dV
class_perf <- mean(Qbar0(class_dV, testdata[, Wnodes]))
joint_dV <- predict(result, newdata = testdata[, Wnodes], pred_fit = "joint")$dV
joint_perf <- mean(Qbar0(joint_dV, testdata[, Wnodes]))
EYd0_perf <- mean(Qbar0(testdata$d0, testdata[, Wnodes]))
c(QaV_perf, class_perf, joint_perf, EYd0_perf)
# perf of true blip approx=0.748

plot(result)

vim <- tx_vim(result)
ggplot(vim, aes(y = node, x = risk_full_fraction, color = model)) + geom_point() + theme_bw() + xlab("VIM")

library(reshape2)
long <- melt(vim, id = c("node", "model"))
ggplot(long, aes(y = node, x = value, color = model)) + geom_point() + facet_wrap(~variable, scales = "free") + 
    theme_bw() 
