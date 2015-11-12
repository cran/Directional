################################
#### BIC to choose the number of components in a 
#### model based clustering using mixtures of von Mises-Fisher distributions
#### Tsagris Michail 4/2015 
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################

bic.mixvmf <- function(x, A = 3) {
  ## x contains the data
  ## A is the maximum number of clusters, set to 3 by default
  x <- as.matrix(x)  ## make sure the x is a matrix
  n <- nrow(x)  ## sample size of the data
  p <- ncol(x)  ## dimensionality of the data
  x <- x/sqrt(rowSums(x^2))  ## makes sure the data are directional
  BIC <- 1:A
  mod <- vmf(x) 
  BIC[1] <- -2 * mod$loglik+ p * log(n)  ## BIC assuming one cluster
  for (vim in 2:A) {
    a <- mix.vmf(x, vim)  ## model based custering for a range of possible clusters
    BIC[vim] <- -2 * a$loglik + ((vim - 1) + vim * p) * log(n)
  }  ## BIC for a range of different clusters
  names(BIC) <- 1:A
  ina <- rep(1, A)
  ina[which.min(BIC)] <- 2  ## chosen number of clusters will 
  ## appear with red on the plot
  plot(1:A, BIC, pch = 10, col = ina, xlab = "Number of components", 
  ylab = "BIC values", type = "b")
  BIC
}