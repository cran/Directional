## Goodness of fit test for grouped data
## The hypothesis is that they come from a
## von Mises-Fisher or a uniform distribution
## May 2016
## References: Arthur Pewsey, Markus Neuhauser, and Graeme D. Ruxton (2013)
## Circular Statistics in R
group.gof <- function(g, ni, m, k, dist = "vm", rads = FALSE, R = 999, ncores = 1) {

  if ( !rads )  m <- m * pi / 180
  d <- length(ni) ##  number of groups
  p <- numeric(d)
  n <- sum(ni)

  if ( dist == "vm" ) {
    for ( i in 1:d )  p[i] <- Directional::pvm( g[i + 1], m, k, rads = TRUE ) -  Directional::pvm( g[i], m, k, rads = TRUE )
  } else if ( dist == "uniform" ) {
    for ( i in 1:d )  p[i] <- ( g[i + 1] - g[i] ) / (2 * pi)
  }

  sj <- cumsum(ni - n * p)
  sjbar <- sum(p * sj)
  ug <- sum(p * (sj - sjbar)^2 ) / n

  if ( ncores <= 1 ) {
    tic <- proc.time()
    ub <- numeric(R)
    for (i in 1:R) {
      x <- Directional::rvonmises(n, m, k, rads = TRUE)
      sim <- as.vector( table( cut(x, breaks = g) ) )
      sjb <- cumsum(sim - n * p)
      sjbarb <- sum(p * sjb)
      ub[i] <- sum( p * (sjb - sjbarb)^2 ) / n
    }
    runtime <- proc.time() - tic

  } else {
    tic <- proc.time()
     cl <- parallel::makePSOCKcluster(ncores)
     doParallel::registerDoParallel(cl)
     ub <- foreach::foreach(vim = 1:R, .combine = rbind, .export = "rvonmises") %dopar% {
       x <- Directional::rvonmises(n, m, k, rads = TRUE)
       sim <- as.vector( table( cut(x, breaks = g) ) )
       sjb <- cumsum(sim - n * p)
       sjbarb <- sum(p * sjb)
       a <- sum( p * (sjb - sjbarb)^2 ) / n
       return(a)
     }
     parallel::stopCluster(cl)
     runtime <- proc.time() - tic
  }

  p.value <- ( sum(ub > ug) + 1 ) / ( R + 1 )
  parameter <- "NA"     ;   names(parameter) <- "df"
  statistic <- ug  ;   names(statistic) <- "Test statistic"
  alternative <- "The distribution is not von Mises-Fisher or uniform distribution"
  method <- "Goodness of fit test for grouped data for testing whether the distribution is von Mises-Fisher or uniform"
  data.name <- c("data")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
}
