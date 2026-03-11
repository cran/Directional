gcpc.means.test <- function(u1, u2, rads = FALSE) {

  lik <- function(para, u1, u2, n1, n2) {
    omega <- para[1] %% (2 * pi)
    r1 <- exp(para[2])  ;  r2 <- exp(para[3])
    g1 <- exp(para[4])  ;  g2 <- exp(para[5])
    phi1 <- u1 - omega  ;  phi2 <- u2 - omega
    a1 <- g1 * cos(phi1)  ;  a2 <- g2 * cos(phi2)
    b1 <- cos(phi1)^2 + sin(phi1)^2/r1
    b2 <- cos(phi2)^2 + sin(phi2)^2/r2
    0.5 * n1 * log(r1) + 0.5 * n2 * log(r2) +
      sum( log( b1 * sqrt(g1^2 + 1) - a1 * sqrt(b1) ) ) +
      sum( log( b2 * sqrt(g2^2 + 1) - a2 * sqrt(b2) ) )
  }

  if ( !rads )  {
    u1 <- u1 * pi/180
    u2 <- u2 * pi/180
  }

  lika <- gcpc.mle(u1, rads = TRUE)
  likb <- gcpc.mle(u2, rads = TRUE)
  lik1 <- lika$loglik + likb$loglik
  g1 <- lika$gamma  ;  g2 <- likb$gamma
  r1 <- lika$rho  ;  r2 <- likb$rho
  m <- 0.5 * (lika$circmu + likb$circmu)
  n1 <- length(u1)  ;  n2 <- length(u2)  ;  n <- n1 + n2
  ini <- c( m, log(r1), log(r2), log(g1), log(g2) )
  lik0 <-  - optim(ini, lik, u1 = u1, u2 = u2, n1 = n1, n2 = n2)$value - n * log(2 * pi)
  stat <- 2 * (lik1 - lik0)
  pvalue <- pchisq(stat, 1, lower.tail = FALSE)

  statistic <- stat  ;   names(statistic) <- "Test statistic"
  parameter <- 1     ;   names(parameter) <- "df"
  alternative <- "The 2 location parameters differ"
  method <- "Location parameter testing using the GCPC distribution"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = pvalue,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)

}
