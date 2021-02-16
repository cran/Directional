################################
#### Summary statistics for cicular data
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
#### References: S Rao Jammalamadaka and A SenGupta (2001)
#### Topics in circular statistics
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 124
################################
circ.summary <- function(u, rads = FALSE, fast = FALSE, tol = 1e-07, plot = TRUE) {
  ## u is an angular variable
  if ( !rads )   u <- u * pi/180
  if (fast) {
    mod <- Rfast::vm.mle(u, tol = tol)
    res <- list(mesos = mod$param[1], kappa = mod$param[2], loglik = mod$loglik)
  } else {
    n <- length(u)  ## sample size
    ## if the data are in degrees we transform them into radians
    ## mesos contains the sample mean
    ## direction
    C <- sum( cos(u) ) / n
    S <- sum( sin(u) )/ n
    Rbar <- sqrt( C^2 + S^2 )  ## mean resultant length
    if (C > 0) {
      mesos <- atan(S/C)
    } else  mesos <- atan(S/C) + pi
    MRL <- Rbar  ## mean resultant length
    circv <- 1 - Rbar
    circs <- sqrt( -2 * log(Rbar) )  ## sample cicrular standard deviation
    ## lik is the von Mises likelihood
    lik <- function(k)   k * sum( cos(u - mesos) ) - n * ( log(besselI( k, 0, expon.scaled = TRUE) ) + k )
	  mod <- optimize(lik, c(0, 50000), maximum = TRUE, tol = tol)
    kappa <- mod$maximum
    ## kappa is the estimated concentration (kappa)
    R <- n * Rbar
    if (Rbar < 2/3) {
      fact <- sqrt(2 * n * ( 2 * R^2 - n * qchisq(0.95, 1) )/ ( R^2 * ( 4 * n - qchisq(0.95, 1)) ) )
      ci <- c(mesos - acos(fact), mesos + acos(fact))
    } else  {
      fact <- sqrt( n^2 - (n^2 - R^2) * exp( qchisq(0.95, 1)/n ) )/R
      ci <- c(mesos - acos(fact), mesos + acos(fact))
    }
    if ( !rads ) {
      mesos <- mesos * 180/pi
      ci <- ci * 180/pi
    }
    if ( plot ) {
      r <- seq(0, 2 * pi, by = 0.01)
      plot(cos(r), sin(r), type = "l", xlab = "Cosinus", ylab = "Sinus", cex.lab = 1.2)
      xx <- seq(-1, 1, by = 0.1)
      yy <- seq(-1, 1, by = 0.1)
      ta <- numeric(length(xx))
      lines(ta, xx, type = "l", lty = 2)
      lines(yy, ta, lty = 2)
      points(cos(u), sin(u))
    }
    res <- list( mesos = mesos, confint = ci, kappa = kappa, MRL = MRL, circvariance = circv, circstd = circs, loglik = mod$objective - n * log(2 * pi) )
  }
  res
}
