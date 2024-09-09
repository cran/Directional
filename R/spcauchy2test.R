spcauchy2test <- function(y1, y2, B = 1) {
  lik0 <- .sp2(y1, y2)$loglik
  lik1 <- Directional::spcauchy.mle(y1)$loglik + Directional::spcauchy.mle(y2)$loglik
  stat <- 2 * lik1 - 2 * lik0
  if ( B == 1 ) {
    pvalue <- pchisq(stat, dim(y1)[2] - 1, lower.tail = FALSE)
  } else {
    n1 <- dim(y1)[1]  ;  n2 <- dim(y2)[1]
    bstat <- numeric(B)
    for (i in 1:B) {
      y1b <- y1[sample(n1, n1, replace = TRUE), ]
      y2b <- y2[sample(n2, n2, replace = TRUE), ]
      lik0 <- .sp2(y1b, y2b)$loglik
      lik1 <- Directional::spcauchy.mle(y1b)$loglik + Directional::spcauchy.mle(y2b)$loglik
      bstat[i] <- lik1 - lik0
    }
    pvalue <- ( sum(bstat >= 0.5 * stat) + 1 ) / (B + 1)
  }

  statistic <- stat  ;   names(statistic) <- "Test statistic"
  parameter <- dim(y1)[2] - 1    ;   names(parameter) <- "df"
  alternative <- "The 2 location parameters differ"
  method <- "Location parameter testing using the spherical Cauchy distribution"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = pvalue,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}




.sp2 <- function(y1, y2, tol = 1e-6) {

  dm <- dim(y1)
  n1 <- dm[1]  ;  d <- dm[2] - 1
  n2 <- dim(y2)[1]

  sp1 <- function(r1, r2, mu, y1, y2, n1, n2, d) {
    a1 <- as.vector(y1 %*% mu)
    a2 <- as.vector(y2 %*% mu)
    n1 * d * log(1 - r1^2) - d * sum( log1p( r1^2 - 2 * r1 * a1 ) ) +
      n2 * d * log(1 - r2^2) - d * sum( log1p( r2^2 - 2 * r2 * a2 ) )
  }

  sp2 <- function(r2, r1, mu, y1, y2, n1, n2, d) {
    a1 <- as.vector(y1 %*% mu)
    a2 <- as.vector(y2 %*% mu)
    n1 * d * log(1 - r1^2) - d * sum( log1p( r1^2 - 2 * r1 * a1 ) ) +
      n2 * d * log(1 - r2^2) - d * sum( log1p( r2^2 - 2 * r2 * a2 ) )
  }

  mu <- Rfast::colmeans( rbind(y1, y2) )
  mu <- mu / sqrt(sum(mu^2) )
  mod1 <- optimize( sp1, c(0, 1), r2 = 0.5, mu = mu, y1 = y1, y2 = y2, n1 = n1, n2 = n2, d = d,
                    maximum = TRUE, tol = 1e-6 )
  mod2 <- optimize( sp2, c(0, 1), r1 = mod1$maximum, mu = mu, y1 = y1, y2 = y2, n1 = n1, n2 = n2, d = d,
                    maximum = TRUE, tol = 1e-6 )
  lik1 <- mod2$objective
  r1 <- mod1$maximum
  r2 <- mod2$maximum

  down <- 1 + r1^2 - 2 * r1 * as.vector( y1 %*% mu)
  mu1 <- Rfast::eachcol.apply(r1 * y1, down, oper = "/")
  down <- 1 + r2^2 - 2 * r2 * as.vector( y2 %*% mu)
  mu2 <- Rfast::eachcol.apply(r2 * y2, down, oper = "/")
  mu <- mu1 + mu2
  mu <- mu / sqrt( sum(mu^2) )

  mod1 <- optimize( sp1, c(0, 1), r2 = r2, mu = mu, y1 = y1, y2 = y2, n1 = n1, n2 = n2, d = d,
                    maximum = TRUE, tol = 1e-6 )
  mod2 <- optimize( sp2, c(0, 1), r1 = mod1$maximum, mu = mu, y1 = y1, y2 = y2, n1 = n1, n2 = n2, d = d,
                    maximum = TRUE, tol = 1e-6 )
  lik2 <- mod2$objective
  r1 <- mod1$maximum
  r2 <- mod2$maximum

  while ( abs( lik2 - lik1 ) > tol ) {
    lik1 <- lik2
    down <- 1 + r1^2 - 2 * r1 * as.vector( y1 %*% mu)
    mu1 <- Rfast::eachcol.apply(r1 * y1, down, oper = "/")
    down <- 1 + r2^2 - 2 * r2 * as.vector( y2 %*% mu)
    mu2 <- Rfast::eachcol.apply(r2 * y2, down, oper = "/")
    mu <- mu1 + mu2
    mu <- mu / sqrt( sum(mu^2) )

    mod1 <- optimize( sp1, c(0, 1), r2 = r2, mu = mu, y1 = y1, y2 = y2, n1 = n1, n2 = n2, d = d,
                      maximum = TRUE, tol = 1e-6 )
    mod2 <- optimize( sp2, c(0, 1), r1 = mod1$maximum, mu = mu, y1 = y1, y2 = y2, n1 = n1, n2 = n2, d = d,
                      maximum = TRUE, tol = 1e-6 )
    lik2 <- mod2$objective
    r1 <- mod1$maximum
    r2 <- mod2$maximum
  }

  n <- n1 + n2
  list(mu = mu, r1 = r1, r2 = r2, loglik = lik2 + n * lgamma( 0.5 * (d + 1) ) - 0.5 * n * (d + 1) * log(pi) - n * log(2) )
}


