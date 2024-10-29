################################
#### Model based clustering using mixtures of von Mises-Fisher distributions
#### Tsagris Michail 4/2015
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################
mixvmf.mle <- function(x, g, n.start = 5, tol = 1e-6, maxiters = 100) {
  ## x contains the data
  ## g is the number of clusters
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  wij <- matrix(nrow = n, ncol = g)
  ka <- numeric(g)
  Apk <- function(p, k)   besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)
  runtime <- proc.time()
  ## Step 1
  l <- 1

  ini <- kmeans(x, g, nstart = n.start)  ## initially a k-means for starting values
  cl <- ini$cluster
  wij <- tabulate(cl)

  if ( min(wij) <= 3 ) {
    mess <- paste( "Too many clusters to fit for this data. Try one less" )
    res <- list(mess = mess, loglik = NA)

  } else {

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
    wlika <- exp(lika)
    rswlika <- Rfast::rowsums(wlika)
    lik[1] <- sum( log( rswlika ) )  ## initial log-likelihood

    ## Step 2
    wij <- wlika / rswlika  ## weights at step 2
    #pj <- Rfast::colmeans(wij) ## weights for step 2

    for (j in 1:g) {
      m1 <- Rfast::eachcol.apply(x, wij[, j])   ## Rfast::colsums(wij[, j] * x)
      mat[j, ] <- m1 / sqrt( sum(m1^2) )  ## mean directions at step 2
      R <- sqrt( sum(m1^2) ) / sum( wij[, j] )  ## mean resultant lengths at step 2
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

    wexplika <- wij * exp(lika)
    rswexplika <- Rfast::rowsums( wexplika )
    lik[2] <- sum( log( rswexplika ) )  ## log-likelihood at step 2
    l <- 2

    ## Step 3 and beyond
    while ( abs(lik[l] - lik[l - 1]) > tol & l < maxiters ) {
      l <- l + 1
      wij <- wexplika / rswexplika  ## weights
      #pj <- Rfast::colmeans(wij)

      for (j in 1:g) {
        m1 <- Rfast::eachcol.apply(x, wij[, j])   ## Rfast::colsums(wij[, j] * x)
        mat[j, ] <- m1 / sqrt( sum(m1^2) )  ## mean directions at step l
        R <- sqrt( sum(m1^2) ) / sum(wij[, j] )  ## mean resultant lengths at step l
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
        log(besselI(ka[j], p/2 - 1, expon.scaled = TRUE) ) - ka[j] + ka[j] * ( x %*% mat[j, ] )
      }

      wexplika <- wij * exp(lika)
      rswexplika <- Rfast::rowsums( wexplika )
      lik[l] <- sum( log( rswexplika ) )
    }  ## log-likelihood at step l

    pj <- Rfast::colmeans(wij)
    loglik <- sum( log( Rfast::colsums( pj * t( exp(lika) ) ) ) )
    ta <- Rfast::rowMaxs(wij)  ## estimated cluster of each observation
    param <- cbind( pj, ka, mat )
    runtime <- proc.time() - runtime
    colnames(param) <- c( "probs", "kappa", paste("mu", 1:p, sep = "") )
    rownames(param) <- paste("Cluster", 1:g, sep = " ")
    res <- list(param = param, loglik = loglik, pred = ta, w = wij, iter = l, runtime = runtime)

  }

  res
}
