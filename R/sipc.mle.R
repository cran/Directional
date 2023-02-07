sipc.mle <- function(x, tol = 1e-6) {

   mag <- function(m, x) {
     a <- as.vector( x %*% m )
     rl <- sum(m^2)
     d <- rl + 1 - a^2
     sqd <- sqrt(d)
     up <- log( (rl + 1) * sqd * ( atan2(sqd, -a) - atan2(sqd, a) + pi ) + 2 * a * d )
     down <- log(d^2)
     - sum(up) + sum(down)
   }
  n <- dim(x)[1]
  mod <- Directional::iag.mle(x)
  da <- nlm( mag, mod$mesi[1, ], x = x, iterlim = 10000 )
  lik1 <-  -da$minimum
  da <- optim( da$estimate, mag, x = x, control = list(maxit = 10000) )
  lik2 <-  -da$value
  while (lik2 - lik1 > tol) {
    lik1 <- lik2
    da <- optim( da$par, mag, x = x, control = list(maxit = 10000) )
    lik2 <-  -da$value
  }

  list(mu = da$par, loglik = lik2 - n * log(4 * pi^2) )
}
