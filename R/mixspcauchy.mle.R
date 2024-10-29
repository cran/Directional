mixspcauchy.mle <- function(x, g, n.start = 5, tol = 1e-6, maxiters = 100) {

  fun2 <- function(wlika, rswlika, mu, x, d, g, lika, tol, maxiters) {

    wij <- wlika / rswlika  ## weights
    wij[wij < 1e-5] <- 1e-5
    wij <- wij / Rfast::rowsums(wij)
    pj <- c( Rfast::colmeans(wij) )   # PANOS
    for (j in 1:g) {
      a2 <- .wspcauchy.wmle(x, w = wij[, j], tol = tol, maxiters = maxiters)
      mu[j, ] <- a2$mesos
      g2 <- sum(mu[j, ]^2)
      a <- as.vector(x %*% mu[j, ])
      #lika[, j] <-  - d * wij[, j] * log( ( sqrt(g2 + 1) - a ) )
      lika[, j] <-  - d * log( ( sqrt(g2 + 1) - a ) ) + log(pj[j])  #PANOS
    }

    #wlika <- wij * exp(lika)
    #rswlika <- Rfast::rowsums(wlika)
    #lik <- sum( log( rswlika ) )

    wlika <-  exp(lika) 		#PANOS
    rswlika <- Rfast::rowsums(wlika) #PANOS
    lik <- sum( log( rswlika ) ) 	#PANOS
    wij <- wlika/rswlika

    list(wij = wij, mu = mu, wlika = wlika, rswlika = rswlika, lika = lika, lik = lik)
  }

  dm <- dim(x)
  n <- dm[1]  ;  d <- dm[2] - 1

  lik <- NULL
  lika <- matrix(nrow = n, ncol = g)
  mu <- matrix(nrow = g, ncol = d + 1)

  runtime <- proc.time()
  ## Step 1

  ini <- kmeans(x, g, nstart = n.start)  ## initially a k-means for starting values
  cl <- ini$cluster
  wij <- tabulate(cl)

  if ( min(wij) <= 3 ) {
    mess <- paste( "Too many clusters to fit for this data. Try one less" )
    res <- list(mess = mess, loglik = NA)

  } else {

    for (j in 1:g) {
      mod <- Directional::spcauchy.mle2(x[cl == j, ])
      gam <- 2 * mod$rho / (1 - mod$rho^2)
      mu[j, ] <- gam * mod$mu
      a <- as.vector( x %*% mu[j, ] )
      lika[, j] <-  - d * log( sqrt(gam^2 + 1) - a )
    }

    wlika <- exp(lika)
    rswlika <- Rfast::rowsums(wlika)
    #lik <- sum( log( rswlika ) )  ## initial log-likelihood

    ep <- fun2(wlika, rswlika, mu, x, d, g, lika, tol, maxiters)
    lik[1] <- ep$lik
    ep2 <- fun2(ep$wlika, ep$rswlika, ep$mu, x, d, g, lika, tol, maxiters)
    lik[2] <- ep2$lik

    i <- 2
    while ( lik[i] - lik[i - 1] > tol & i < maxiters ) {
      i <- i + 1
      ep <- ep2
      ep2 <- fun2(ep$wlika, ep$rswlika, ep$mu, x, d, g, lika, tol, maxiters)
      lik[i] <- ep2$lik
    }
    res <- ep2
    if ( ep$lik > ep2$lik )  res <- ep

    pj <- Rfast::colmeans(res$wij)
    loglik <- res$lik   #sum( log( Rfast::colsums( pj * t( exp( res$lika ) ) ) ) )
    ta <- Rfast::rowMaxs(res$wij)  ## estimated cluster of each observation
    gama <- sqrt( Rfast::rowsums(mu^2) )
    dirmu <- mu / gama
    rho <- ( sqrt(gama^2 + 1) - 1 ) / gama
    param <- cbind(pj, gama, mu)
    dirparam <- cbind(pj, rho, dirmu)

    runtime <- proc.time() - runtime
    colnames(param) <- c( "probs", "gama", paste("mu", 1:(d + 1), sep = "") )
    colnames(dirparam) <- c( "probs", "rho", paste("mu", 1:(d + 1), sep = "") )
    rownames(param) <- rownames(dirparam) <- paste("Cluster", 1:g, sep = " ")

    #nz <- res$wij
    #ind <- 1:g
    #thresh <- -744
    #for (i in 1:n) {
    #	index <- ind[nz[i, ] < exp(thresh)]
    #	nz[i, index] <- rep( exp(thresh), length(index) )
    #	nz[i, ] <- nz[i, ]/sum(nz[i, ])
    #}
    #wij <- nz

    res <- list( param = param, dirparam = dirparam,
                 loglik = loglik + n * lgamma( 0.5 * (d + 1) ) - 0.5 * n * (d + 1) * log(pi) - n * log(2),
                 pred = ta, w = res$wij, iters = i, runtime = runtime )
  }

  res

}







.wspcauchy.wmle <- function(x, w, tol = 1e-6, maxiters = 100) {

  wx <- w * x
  sw <- sum(w)
  dm <- dim(x)
  n <- dm[1]  ;  d <- dm[2] - 1

  sp <- function(rho, mu, x, n, d, w, sw) {
    a <- as.vector(x %*% mu)
    d * sw * log(1 - rho^2) - d * sum( w * log1p( rho^2 - 2 * rho * a ) )
  }

  mu <- Rfast::colmeans(wx)
  mu <- mu / sqrt(sum(mu^2) )
  mod <- optimize(sp, c(0, 1), mu = mu, x = x, n = n, d = d, w = w, sw = sw, maximum = TRUE, tol = 1e-6 )
  rho <- mod$maximum
  lik1 <- mod$objective

  down <- 1 + rho^2 - 2 * rho * as.vector( x %*% mu)
  mu <- Rfast::eachcol.apply(rho * wx, down, oper = "/")
  mu <- mu / sqrt( sum(mu^2) )
  mod <- optimize(sp, c(0, 1), mu = mu, x = x, n = n, d = d, w = w, sw = sw, maximum = TRUE, tol = 1e-6 )
  rho <- mod$maximum
  lik2 <- mod$objective
  i <- 2
  while ( abs( lik2 - lik1 ) > tol & i < maxiters ) {
    i <- i + 1
    lik1 <- lik2
    down <- 1 + rho^2 - 2 * rho * as.vector( x %*% mu)
    mu <- Rfast::eachcol.apply(rho * wx, down, oper = "/")
    mu <- mu / sqrt( sum(mu^2) )
    mod <- optimize(sp, c(0, 1), mu = mu, x = x, n = n, d = d, w = w, sw = sw, maximum = TRUE, tol = 1e-6 )
    rho <- mod$maximum
    lik2 <- mod$objective
  }
  gama <- 2 * rho / (1 - rho^2)
  list( mesos = gama * mu, loglik = lik2 )
}


