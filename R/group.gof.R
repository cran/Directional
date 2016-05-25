## Goodness of fit test for grouped data
## The hypothesis is that they come from a
## von Mises-Fisher distribution
## May 2016
## References: Arthur Pewsey, Markus Neuhauser, and Graeme D. Ruxton (2013)
## Circular Statistics in R


group.gof <- function(g, ni, m, k, dist = "vm", rads = FALSE, R = 999, ncores = 1) {

  if ( rads == FALSE )  {
     m <- m * pi / 180
  }

  d <- length(ni) ##  number of groups
  p <- numeric(d)
  n <- sum(ni)

  if ( dist == "vm" ) {
    for ( i in 1:d ) {
      p[i] <- pvm( g[i + 1], m, k, rads = TRUE ) -  pvm( g[i], m, k, rads = TRUE )
    }

  } else if ( dist == "uniform" ) {
      p[i] <- ( g[i + 1] - g[i] ) / (2 * pi)
  }

  sj <- cumsum(ni - n * p)
  sjbar <- sum(p * sj)
  ug <- 1/n * sum(p * (sj - sjbar)^2 )

  if (ncores <= 1) {

    tic <- proc.time()

    ub <- numeric(R)

    for (i in 1:R) {
      x <- rvonmises(n, m, k, rads = TRUE)
      sim <- as.vector( table( cut(x, breaks = g) ) )
      sjb <- cumsum(sim - n * p)
      sjbarb <- sum(p * sjb)
      ub[i] <- 1/n * sum( p * (sjb - sjbarb)^2 )
    }

    runtime <- proc.time() - tic

  } else {

    tic <- proc.time()

     cl <- makePSOCKcluster(ncores)
     registerDoParallel(cl)
     ub <- foreach(vim = 1:R, .combine = rbind, .export = "rvonmises") %dopar% {
       x <- rvonmises(n, m, k, rads = TRUE)
       sim <- as.vector( table( cut(x, breaks = g) ) )
       sjb <- cumsum(sim - n * p)
       sjbarb <- sum(p * sjb)
       a <- 1/n * sum( p * (sjb - sjbarb)^2 )
       return(a)
     }

     stopCluster(cl)

     runtime <- proc.time() - tic

  }

  pval <- ( sum(ub > ug) + 1 ) / ( R + 1 )

  info <- c(ug, pval)
  names(info) <- c("statistic", "p-value")
  list(info = info, runtime = runtime)

}
