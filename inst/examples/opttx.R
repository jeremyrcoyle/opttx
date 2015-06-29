library(opttx)
Qbar0 <- function(A, W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    W4 <- W[, 4]
    Qbar <- ifelse(W4 > 0, plogis(1 - W1^2 + 3 * W2 + A * (5 * W3^2 - 4.45)), plogis(-0.5 - W3 + 2 * W1 * W2 + A * (3 * 
        abs(W2) - 1.5)))
    return(Qbar)
}

g0 <- function(W) {
    W1 <- W[, 1]
    W2 <- W[, 2]
    W3 <- W[, 3]
    W4 <- W[, 4]
    
    # rep(0.5, nrow(W))
    plogis(0.25 * W1 - 0.1 * W2)
}

gen_data <- function(n = 1000, p = 4) {
    W <- matrix(rnorm(n * p), nrow = n)
    colnames(W) <- paste("W", seq_len(p), sep = "")
    A <- rbinom(n, 1, g0(W))
    u <- runif(n)
    Y <- as.numeric(u < Qbar0(A, W))
    Y0 <- as.numeric(u < Qbar0(0, W))
    Y1 <- as.numeric(u < Qbar0(1, W))
    d0 <- as.numeric(Qbar0(1, W) > Qbar0(0, W))
    Yd0 <- as.numeric(u < Qbar0(d0, W))
    data.frame(W, A, Y, Y0, Y1, Yd0, d0, blip = Qbar0(1, W) - Qbar0(0, W))
}

data <- gen_data(1000, 5)

system.time({
    result <- opt_tmle(data, blip_type = "cl.surlog")
})

# perf of true blip approx=0.595
mean(pmax(Qbar0(0, data[, result$nodes$Wnodes]), Qbar0(1, data[, result$nodes$Wnodes])))

# confirm D1 isn't too crazy quantile(extract_vals(result$folds, result$split_preds)$D1)
# quantile(cv_predict_original(result$fits$blip_fit))

print(result)
plot(result)

vim <- tx_vim(result)
ggplot(vim, aes(y = node, x = risk_full_fraction, color = model)) + geom_point() + theme_bw() + xlab("VIM")

library(reshape2)
long <- melt(vim, id = c("node", "model"))
ggplot(long, aes(y = node, x = value, color = model)) + geom_point() + facet_wrap(~variable, scales = "free") + theme_bw() 

