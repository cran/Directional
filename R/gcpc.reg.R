gcpc.reg <- function(y, x, rads = TRUE, reps = 20, xnew = NULL) {

  lik <- function(param, y, x, y1, y2, y12, n, rho) {
    be <- matrix(param, ncol = 2)
    mu <- x %*% be
    g2 <- Rfast::rowsums(mu^2)
    ksi <- mu / sqrt(g2)
    s1 <- ksi[, 1]^2 + ksi[, 2]^2/rho
    s12 <- ksi[, 1] * ksi[, 2] * (1 - 1/rho)
    s2 <- ksi[, 2]^2 + ksi[, 1]^2/rho
    a <- Rfast::rowsums(y * mu)
    B <- y1 * s1 + 2 * y12 * s12 + y2 * s2
    n * 0.5 * log(rho) + sum( log( B * sqrt(g2 + 1) - a * sqrt(B) ) )
  }

  lik1 <- function(rho, mu, g2, ksi, a, y, x, y1, y2, y12, n) {
    rho <- 1 / ( 1 + exp(-rho) )
    s1 <- ksi[, 1]^2 + ksi[, 2]^2/rho
    s12 <- ksi[, 1] * ksi[, 2] * (1 - 1/rho)
    s2 <- ksi[, 2]^2 + ksi[, 1]^2/rho
    B <- y1 * s1 + 2 * y12 * s12 + y2 * s2
    n * 0.5 * log(rho) + sum( log( B * sqrt(g2 + 1) - a * sqrt(B) ) )
  }

  likreg <- function(param, y, x, y1, y2, y12, n) {
    rho <- 1 / ( 1 + exp(-param[1]) )
    be <- matrix(param[-1], ncol = 2)
    mu <- x %*% be
    g2 <- Rfast::rowsums(mu^2)
    ksi <- mu / sqrt(g2)
    s1 <- ksi[, 1]^2 + ksi[, 2]^2/rho
    s12 <- ksi[, 1] * ksi[, 2] * (1 - 1/rho)
    s2 <- ksi[, 2]^2 + ksi[, 1]^2/rho
    a <- Rfast::rowsums(y * mu)
    B <- y1 * s1 + 2 * y12 * s12 + y2 * s2
    n * 0.5 * log(rho) + sum( log( B * sqrt(g2 + 1) - a * sqrt(B) ) )
  }


  if ( !is.matrix(y) ) {
    if ( !rads )  y <- y * pi/180
    y <- cbind( cos(y), sin(y) )
  }

  x <- model.matrix(y~., as.data.frame(x), )
  n <- dim(y)[1]  ;  p <- dim(x)[2]
  runtime <- proc.time()
  y1 <- y[, 1]^2  ;  y2 <- y[, 2]^2  ;  y12 <- y[, 1] * y[, 2]

  be <- as.vector( spml.reg(y, x[, -1], rads = TRUE)$be )
  be <- matrix( optim(be, lik, y = y, x = x, y1 = y1, y2 = y2, y12 = y12, n = n, rho = 0.5)$par, ncol = 2 )
  mu <- x %*% be
  g2 <- Rfast::rowsums(mu^2)
  ksi <- mu / sqrt(g2)
  a <- Rfast::rowsums(y * mu)
  modrho <- optimise( lik1, c(0.001, 1000), mu = mu, g2 = g2, ksi = ksi, a = a, y = y, x = x, y1 = y1,
                      y2 = y2, y12 = y12, n = n, maximum = TRUE )
  rho <- modrho$maximum
  be <- as.vector( optim(as.vector(be), lik, y = y, x = x, y1 = y1, y2 = y2, y12 = y12, n = n, rho = rho)$par )
  mod <- optim( c( log(rho / (1 - rho)), be ), likreg, y = y, x = x, y1 = y1, y2 = y2,
                y12 = y12, n = n, control = list(maxit = 5000) )
  lika <- mod$value
  mod <- optim( mod$par, likreg, y = y, x = x, y1 = y1, y2 = y2, y12 = y12, n = n,
                control = list(maxit = 5000) )
  likb <- mod$value
  while ( lika - likb > 1e-6 ) {
    lika <- likb
    mod <- optim( mod$par, likreg, y = y, x = x, y1 = y1, y2 = y2, y12 = y12, n = n,
                control = list(maxit = 5000), hessian = TRUE )
    likb <- mod$value
  }

  se <- solve(mod$hessian)
  runtime <- proc.time() - runtime
  rho <- 1 / ( 1 + exp(-mod$par[1]) )
  be <- matrix(mod$par[-1], ncol = 2)
  serho <- rho * (1 - rho) * sqrt( se[1, 1] )
  seb <- matrix( sqrt( diag(se)[-1] ), ncol = 2 )

  colnames(be) <- c("Cosinus of y", "Sinus of y")
  rownames(be) <- colnames(x)
  colnames(seb) <- c("Cosinus of y", "Sinus of y")
  rownames(seb) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix( ~., data.frame(xnew) )
    est <- xnew %*% be
    est <- ( atan(est[, 2] / est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
    if ( !rads )  est <- est * 180 / pi
  }

  list( runtime = runtime, beta = be, seb = seb, rho = rho, serho = serho,
        loglik = -mod$value - n * log(2 * pi), est = est )
}

