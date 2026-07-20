################################
#### BIC to choose the number of components in a
#### model based clustering using mixtures of von Mises-Fisher distributions
#### Tsagris Michail 4/2015
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################
bic.mixvmf <- function(x, G = 5, n.start = 5, tol = 1e-6, maxiters = 500) {
  ## x contains the data
  ## A is the maximum number of clusters, set to 3 by default
  runtime <- proc.time()
  logn <- log( dim(x)[1] )  ## sample size of the data
  p <- dim(x)[2]  ## dimensionality of the data
  bic <- icl <- 1:G
  mod <- Directional::vmf.mle(x)
  bic[1] <- icl[1] <-  - 2 * mod$loglik+ p * logn  ## BIC assuming one cluster
  for (vim in 2:G) {
    a <- Directional::mixvmf.mle(x, vim, n.start = n.start, tol = tol, maxiters = maxiters)  ## model based clustering for some possible clusters
    if  ( !is.na(a$loglik) ) {
      bic[vim] <-  - 2 * a$loglik + ( (vim - 1) + vim * p ) * logn
      icl[vim] <- bic[vim] - 2 * sum( a$w * log(a$w), na.rm = TRUE )
    } else  {
      bic[vim] <- bic[vim - 1]
      icl[vim] <- icl[vim - 1]
    }
  }  ## BIC for a range of different clusters
  runtime <- proc.time() - runtime
  names(bic) <- 1:G
  ina <- rep(1, G)
  ina[ which.min(icl) ] <- 3  ## chosen number of clusters will appear with red on the plot
  plot(1:G, icl, col = ina, xlab = "Number of components", ylab = "ICL values", cex.lab = 1.3, cex.axis = 1.3)
  abline(v = 1:G, lty = 2, col = "lightgrey")
  abline(h = seq(min(icl, na.rm = FALSE), max(icl, na.rm = FALSE), length = 10), lty = 2, col = "lightgrey" )
  lines(1:G, icl, lwd = 2)
  points(1:G, icl, pch = 9, col = ina)
  list(bic = bic, icl = icl, runtime = runtime)
}
