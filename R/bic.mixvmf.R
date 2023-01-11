################################
#### BIC to choose the number of components in a
#### model based clustering using mixtures of von Mises-Fisher distributions
#### Tsagris Michail 4/2015
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################
bic.mixvmf <- function(x, G = 5, n.start = 20) {
  ## x contains the data
  ## A is the maximum number of clusters, set to 3 by default
  runtime <- proc.time()
  logn <- log( dim(x)[1] )  ## sample size of the data
  p <- dim(x)[2]  ## dimensionality of the data
  bic <- 1:G
  mod <- Directional::vmf.mle(x)
  bic[1] <-  - 2 * mod$loglik+ p * logn  ## BIC assuming one cluster
  for (vim in 2:G) {
    a <- Directional::mixvmf.mle(x, vim, n.start = n.start)  ## model based clustering for some possible clusters
    bic[vim] <-  -2 * a$loglik + ( (vim - 1) + vim * p ) * logn
  }  ## BIC for a range of different clusters
  runtime <- proc.time() - runtime
  names(bic) <- 1:G
  ina <- rep(1, G)
  ina[which.min(bic)] <- 3  ## chosen number of clusters will
  ## appear with red on the plot
  plot(1:G, bic, col = ina, xlab = "Number of components",
  ylab = "BIC values")
  abline(v = 1:G, lty = 2, col = "lightgrey")
  abline(h = seq(min(bic), max(bic), length = 10), lty = 2, col = "lightgrey" )
  lines(1:G, bic, lwd = 2)
  points(1:G, bic, pch = 9, col = ina)
  list(bic = bic, runtime = runtime)
}
