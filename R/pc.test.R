pc.test <- function(x, B = 1, tol = 1e-6) {
  mod <- Directional::sespc.mle(x, tol = tol)
  stat <- 2 * mod$loglik - 2 * mod$sipc.loglik
  if (B == 1) {
    pvalue <- pchisq(stat, 2, lower.tail = FALSE)
    res <- c(stat, pvalue)
    names(res) <- c("test", "p-value")
  } else {
    tb <- numeric(B)
    n <- dim(x)[1]
    for (i in 1:B) {
      nu <- Rfast2::Sample.int(n, n, replace = TRUE)
      mod <- Directional::sespc.mle(x[nu, ], tol = tol)
      tb[i] <- 2 * mod$loglik - 2 * mod$sipc.loglik
    }
    res <- c( stat, (sum(tb > stat) + 1) / (B + 1) )
    names(res) <- c('test', 'Bootstrap p-value')
  }
  res
}
