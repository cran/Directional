bic.mixpkbd <- function(x, G = 5, n.start = 5, tol = 1e-6, maxiters = 500) {
    ## x contains the data
    ## A is the maximum number of clusters, set to 3 by default
    runtime <- proc.time()
    logn <- log( dim(x)[1] )  ## sample size of the data
    p <- dim(x)[2]  ## dimensionality of the data
    bic <- icl <- 1:G
    mod <- Directional::pkbd.mle(x)
    bic[1] <- icl[1] <-  - 2 * mod$loglik+ p * logn  ## BIC assuming one cluster
    for (vim in 2:G) {
      a <- Directional::mixpkbd.mle(x, vim, n.start = n.start, tol = tol, maxiters = maxiters)  ## model based clustering for some possible clusters
      bic[vim] <-  -2 * a$loglik + ( (vim - 1) + vim * p ) * logn
      icl[vim] <- bic[vim] - sum( a$w * log(a$w), na.rm = TRUE )
    }  ## BIC for a range of different clusters
    runtime <- proc.time() - runtime
    names(bic) <- 1:G
    ina <- rep(1, G)
    ina[ which.min(bic) ] <- 3  ## chosen number of clusters will appear with red on the plot
    plot(1:G, bic, col = ina, xlab = "Number of components", ylab = "BIC values", cex.lab = 1.3, cex.axis = 1.3)
    abline(v = 1:G, lty = 2, col = "lightgrey")
    abline(h = seq(min(bic, na.rm = FALSE), max(bic, na.rm = FALSE), length = 10), lty = 2, col = "lightgrey" )
    lines(1:G, bic, lwd = 2)
    points(1:G, bic, pch = 9, col = ina)
    list(bic = bic, icl = icl, runtime = runtime)
  }
