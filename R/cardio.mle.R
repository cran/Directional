cardio.mle <- function(x, rads = FALSE) {
  if ( !rads )   x <- x * pi/180
  cardio <- function(para, x) {
    rho2 <- cos(para[1])
    mu <- para[2] %% pi
    - sum( log1p( rho2 * cos(x - mu) ) )
  }
  mod <- optim( c( runif(1), mean(x) ), cardio, x = x, control = list(maxit = 5000) )
  rho <- 0.5 * cos( mod$par[1] )
  mu <- mod$par[2] %% pi
  param <- c(mu, rho)
  names(param) <- c("mu", "rho")
  n <- length(x)
  list( loglik = n * log(2 * pi) - mod$value, param = param)
}
