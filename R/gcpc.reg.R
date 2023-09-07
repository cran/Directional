gcpc.reg <- function(y, x, rads = TRUE, reps = 20, xnew = NULL) {

  lik <- function(param, y, x, y1, y2, y12, n, rho) {
    be <- matrix(param, ncol = 2)
    mu <- x %*% be
    g2 <- rowsums(mu^2)
    ksi <- mu / sqrt(g2)
    s1 <- ksi[, 1]^2 + ksi[, 2]^2/rho
    s12 <- ksi[, 1] * ksi[, 2] * (1 - 1/rho)
    s2 <- ksi[, 2]^2 + ksi[, 1]^2/rho
    a <- rowsums(y * mu)
    B <- y1 * s1 + 2 * y12 * s12 + y2 * s2
    n * 0.5 * log(rho) + sum( log( B * sqrt(g2 + 1) - a * sqrt(B) ) )
  }

  lik1 <- function(param, y, x, y1, y2, y12, n) {
    rho <- 1 / ( 1 + exp(-param[1]) )
    be <- matrix(param[-1], ncol = 2)
    mu <- x %*% be
    g2 <- rowsums(mu^2)
    ksi <- mu / sqrt(g2)
    s1 <- ksi[, 1]^2 + ksi[, 2]^2/rho
    s12 <- ksi[, 1] * ksi[, 2] * (1 - 1/rho)
    s2 <- ksi[, 2]^2 + ksi[, 1]^2/rho
    a <- rowsums(y * mu)
    B <- y1 * s1 + 2 * y12 * s12 + y2 * s2
    n * 0.5 * log(rho) + sum( log( B * sqrt(g2 + 1) - a * sqrt(B) ) )
  }

  lik2 <- function(rho, reps = reps) {
    a <- matrix(nrow = reps, ncol = 2 * p + 1)
    for (i in 1:reps) {
      mod <- optim( rnorm(2 * p), lik, y = y, x = x, y1 = y1, y2 = y2,
                    y12 = y12, n = n, rho = rho, control = list(maxit = 5000) )
      lika1 <-  -mod$value
      mod <- optim( mod$par, lik, y = y, x = x, y1 = y1, y2 = y2,
                    y12 = y12, n = n, rho = rho, control = list(maxit = 5000) )
      lika2 <-  -mod$value
      while (lika2 - lika1 > 1e-6) {
        lika1 <- lika2
        mod <- optim( mod$par, lik, y = y, x = x, y1 = y1, y2 = y2,
                      y12 = y12, n = n, rho = rho, control = list(maxit = 5000) )
        lika2 <- -mod$value
      }
      mod <-  optim( mod$par, lik, y = y, x = x, y1 = y1, y2 = y2,
                     y12 = y12, n = n, rho = rho, control = list(maxit = 5000) )
      a[i, ] <- c( -mod$value, mod$par )
    }  ## end  for (i in 1:reps) {
    ind <- which.max(a[, 1])
    a[ind, ]
  }

  if ( !is.matrix(y) ) {
    if ( !rads )  y <- y * pi/180
    y <- cbind( cos(y), sin(y) )
  }

  x <- model.matrix(y~., as.data.frame(x), )
  n <- dim(y)[1]  ;  p <- dim(x)[2]
  runtime <- proc.time()
  y1 <- y[, 1]^2  ;  y2 <- y[, 2]^2  ;  y12 <- y[, 1] * y[, 2]

  likreg <- function(rho)  lik2(rho, reps = reps)[1]
  rho <- optimise(likreg, c(0.001, 1), maximum = TRUE)$maximum
  par <- lik2(rho, reps = reps)[-1]
  par <- c( log( rho/(1 - rho) ), par)
  mod <- optim( par, lik1, y = y, x = x, y1 = y1, y2 = y2,
                y12 = y12, n = n, control = list(maxit = 5000), hessian = TRUE )
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

