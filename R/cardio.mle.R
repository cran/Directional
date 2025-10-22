cardio.mle <- function(x, rads = FALSE) {
  if ( !rads )   x <- x * pi/180
  n <- length(x)

  cardio <- function(para, x) {
    rho <- para[1]
    mu <- para[2]
    - sum( log1p( 2 * rho * cos(x - mu) ) )
  }

  C <- sum( cos(x) ) / n  ;  S <- sum( sin(x) )/ n
  rho <- sqrt( C^2 + S^2 )  ## mean resultant length
  if (rho > 0.5)  rho <- 0.5
  mu <- ( atan(S / C) + pi * I(C < 0) ) %% (2 * pi)
  mod <- optim( c( rho, mu), cardio, x = x, control = list(maxit = 5000),
                method = "L-BFGS-B", lower = c(0, 0), upper = c(0.5, 2 * pi)  )
  rho <- mod$par[1]
  mu <- mod$par[2]
  param <- c(mu, rho)
  names(param) <- c("mu", "rho")
  list( loglik = -n * log(2 * pi) - mod$value, param = param)
}
