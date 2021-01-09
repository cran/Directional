################################
#### Model based clustering using mixtures of von Mises-Fisher distributions
#### Tsagris Michail 4/2015
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################
mixvmf.mle <- function(x, g, n.start = 10) {
  ## x contains the data
  ## g is the number of clusters
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  pij <- matrix(nrow = n, ncol = g)
  ka <- numeric(g)
  Apk <- function(p, k)   besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)
  runtime <- proc.time()
  ## Step 1
  l <- 1

  ini <- kmeans(x, g, nstart = n.start)  ## initially a k-means for starting values
  cl <- ini$cluster
  w <- tabulate(cl)

  if ( min(w) <= 3 ) {
    mess <- paste( "Too many clusters to fit for this data. Try one less" )
    res <- list(mess = mess, loglik = NA)

  } else {
    w <- w/n  #'# initial weights
    m1 <- ini$centers
    Rk <- sqrt( Rfast::rowsums(m1^2) )  ## mean resultant lengths of the initical clusters
    mat <- m1/Rk  ## initial mean directions

    for (j in 1:g) {
      R <- Rk[j]
      k1 <- R * (p - R^2)/(1 - R^2)
      apk <- Apk(p, k1)
      k2 <- k1 - ( apk - R)/( 1 - apk^2 - (p - 1)/k1 * apk )
      while (abs(k2 - k1) > 1e-07) {
        k1 <- k2
        apk <- Apk(p, k1)
        k2 <- k1 - (apk - R)/( 1 - apk^2 - (p - 1)/k1 * apk )
      }
      ka[j] <- k2  ## initial concentration parameters
      lika[, j] <- (p/2 - 1) * log(ka[j]) - 0.5 * p * log(2 * pi) -
      log( besselI(ka[j], p/2 - 1, expon.scaled = TRUE) ) - ka[j] + ka[j] * (x %*% mat[j, ])
    }
    wlika <- w * exp(lika)
    rswlika <- Rfast::rowsums(wlika)
    lik[1] <- sum( log( rswlika ) )  ## initial log-likelihood

    l <- 2
    ## Step 2
    pij <- wlika / rswlika  ## weights at step 2
    w <- Rfast::colmeans(pij) ## weights for step 2

    for (j in 1:g) {
      m1 <- Rfast::eachcol.apply(x, pij[, j])   ## Rfast::colsums(pij[, j] * x)
      mat[j, ] <- m1 / sqrt( sum(m1^2) )  ## mean directions at step 2
      R <- sqrt( sum(m1^2) ) / sum( pij[, j] )  ## mean resultant lengths at step 2
      k1 <- R * (p - R^2)/(1 - R^2)
      apk <- Apk(p, k1)
      k2 <- k1 - ( apk - R)/( 1 - apk^2 - (p - 1)/k1 * apk )
      while (abs(k2 - k1) > 1e-07) {
        k1 <- k2
        apk <- Apk(p, k1)
        k2 <- k1 - (apk - R)/( 1 - apk^2 - (p - 1)/k1 * apk )
      }
      ka[j] <- k2
      lika[, j] <- (p/2 - 1) * log(ka[j]) - 0.5 * p * log(2 * pi) -
      log(besselI(ka[j], p/2 - 1, expon.scaled = TRUE) ) - ka[j] + ka[j] * (x %*% mat[j, ])
    }

    wexplika <- w * exp(lika)
    lik[2] <- sum( log( Rfast::rowsums( wexplika ) ) )  ## log-likelihood at step 2
    ## Step 3 and beyond
    while ( lik[l] - lik[l - 1] > 1e-05 ) {
      l <- l + 1
      pij <- wexplika / Rfast::rowsums( wexplika )  ## weights
      w <- Rfast::colmeans(pij)

      for (j in 1:g) {
        m1 <- Rfast::eachcol.apply(x, pij[, j])   ## Rfast::colsums(pij[, j] * x)
        mat[j, ] <- m1 / sqrt( sum(m1^2) )  ## mean directions at step l
        R <- sqrt( sum(m1^2) ) / sum(pij[, j])  ## mean resultant lengths at step l
        k1 <- R * (p - R^2)/(1 - R^2)
        apk <- Apk(p, k1)
        k2 <- k1 - (apk - R)/( 1 - apk^2 - (p - 1)/k1 * apk )
        while (abs(k2 - k1) > 1e-07) {
          k1 <- k2
          apk <- Apk(p, k1)
          k2 <- k1 - (apk - R)/( 1 - apk^2 - (p - 1)/k1 * apk )
        }
        ka[j] <- k2
        lika[, j] <- (p/2 - 1) * log(ka[j]) - 0.5 * p * log(2 * pi) -
        log(besselI(ka[j], p/2 - 1, expon.scaled = TRUE) ) - ka[j] + ka[j] * (x %*% mat[j, ])
      }

      wexplika <- w * exp(lika)
      lik[l] <- sum( log( Rfast::rowsums( wexplika ) ) )
    }  ## log-likelihood at step l

    ta <- Rfast::rowMaxs(pij)  ## estimated cluster of each observation
    param <- cbind( mat, ka, table(ta)/n )
    runtime <- proc.time() - runtime
    colnames(param) <- c( paste("mu", 1:p, sep = ""), 'kappa', 'probs' )
    rownames(param) <- paste("Cluster", 1:g, sep = " ")
    res <- list(param = param, loglik = lik[l], pred = ta, iter = l, runtime = runtime)

  }
  res
}
