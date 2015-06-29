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
    rep(0.5, nrow(W))
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
    result <- opt_tmle(data)
})

print(result)
plot(result) 
