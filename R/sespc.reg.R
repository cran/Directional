sespc.reg <- function(y, x, con = TRUE, xnew = NULL, lati = 10, longi = 10, tol = 1e-06) {

  lati <- seq(0, 360, len = lati) * pi / 180
  longi <- seq(0, 360, len = longi) * pi / 180
  x1 <- cos(lati) * cos(longi)
  x2 <- cos(lati) * sin(longi)
  x3 <- sin(lati)
  lat <- asin( x3 / sqrt(x1^2 + x2^2 + x3^2) ) * 180 / pi
  long <- atan2(x2, x1) * 180 / pi
  wa <- Directional::euclid( expand.grid(lat, long) )
  wa <- round(wa, 10)
  wa <- unique(wa)

  my <- Directional::sespc.mle(y)$mu
  my <- my / sqrt( sum(my^2) )

  B <- dim(wa)[1]
  mod <- list()
  lik <- numeric(B)
  x <- model.matrix( ~., data.frame(x) )

  for (i in 1:B) {
    rot <- t( Directional::rotation(my, wa[i, ]) )
    y1 <- y %*% rot
    mod[[ i ]] <- .sespc.reg2(y1, x, con = con, xnew = xnew, tol = tol)
    lik[i] <- mod[[ i ]]$loglik
  }

  ep <- which.max(lik)
  mod[[ ep ]]

}


.sespc.reg2 <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-06) {
  ## y is the spherical data,  unit vectors
  ## x is the independent variable(s)
  ## if you want a constant term, set con eaual to TRUE (default value)
  if ( !con ) x <- x[, -1]
  n <- dim(y)[1]
  y1 <- y[, 1]^2   ;  y2 <- y[, 2]^2   ;  y3 <- y[, 3]^2
  y12 <- y[, 1] * y[, 2]  ;    y13 <- y[, 1] * y[, 3]
  y23 <- y[, 2] * y[, 3]
  za <- cbind(y1, y12, y13, y12, y2, y23, y13, y23, y3)
  ### inner function for SESPC regression
  reg2 <- function(param, y, x, za) {
    the1 <- param[1]
    the2 <- param[2]
    heta <- sqrt( the1^2 + the2^2 + 1 ) - 1
    be <- matrix(param[ - c( 1:2 ) ], ncol = 3)
    m <- x %*% be
    m0 <- sqrt( m[, 2]^2 + m[, 3]^2 )
    rl <- Rfast::rowsums(m^2)
    x1b <- cbind( -m0^2, m[, 1] * m[, 2], m[, 1] * m[, 3] ) / ( m0 * sqrt(rl) )
    x2b <- cbind( 0, -m[, 3], m[, 2] ) / m0
    z12 <- x1b[, 1] * x1b[, 2]  ;    z13 <- x1b[, 1] * x1b[, 3]
    z23 <- x1b[, 2] * x1b[, 3]
    tx1 <- cbind(x1b[, 1]^2, z12, z13, z12, x1b[, 2]^2, z23, z13, z23, x1b[, 3]^2)
    z23 <- x2b[, 2] * x2b[, 3]
    tx2 <- cbind(0, 0, 0, 0, x2b[, 2]^2, z23, 0, z23, x2b[, 3]^2)
    vinv <- cbind(x1b[, 1] * x2b, x1b[, 1] * x2b[, 2], 2 * x1b[, 2] * x2b[, 2], x1b[, 3] * x2b[, 2] + x1b[, 2] * x2b[, 3],
                  x1b[, 1] * x2b[, 3], x1b[, 3] * x2b[, 2] + x1b[, 2] * x2b[, 3], 2 * x1b[, 3] * x2b[, 3]  )
    a <- Rfast::rowsums( y * m )
    b <- 1 + Rfast::rowsums( (the2 * vinv + the1 * (tx1 - tx2) + heta * (tx1 + tx2) ) * za )
    E <- b * rl + b - a^2
    sqe <- sqrt(E)
    up <- log( b * (rl + 1) * sqe * ( atan2(sqe, -a) - atan2(sqe, a) + pi ) + 2 * a * E )
    down <- log( b * E^2 )
    - sum(up) + sum(down)
  }

  ini <- rnorm( 3 * dim(x)[2] + 2 )
  suppressWarnings({
    val1 <- nlm(reg2, ini, y = y, x = x, za = za, iterlim = 5000)
    val2 <- nlm(reg2, val1$estimate, y = y, x = x, za = za, iterlim = 5000)
    while (val1$minimum - val2$minimum > 1e-06) {
      val1 <- val2
      val2 <-  nlm(reg2, val1$estimate, y = y, x = x, za = za, iterlim = 5000)
    }
    da <- optim(val2$estimate, reg2, y = y, x = x, za = za, control = list(maxit = 10000), hessian = TRUE)
  })
  the1 <- da$par[1]
  the2 <- da$par[2]
  rho <- sqrt( the1^2 + the2^2 + 1 ) - sqrt( the1^2 + the2^2 )
  psi <- 0.5 * acos( 2 * the1 / ( 1/rho - rho ) )
  be <- matrix(da$par[ -c( 1:2 ) ], ncol = 3)
  m <- x %*% be
  rl <- sqrt( Rfast::rowsums(m^2) )
  est <- m/rl
  #di <- sum( y * est )
  di <- Rfast::XopY.sum(y, est) 
  se <- sqrt( diag( solve(da$hessian) ) )
  seb <- matrix(se[-c(1:2)], ncol = 3)

  if ( is.null(xnew) ) {
    est <- est
  } else {
    xnew <- model.matrix( ~., data.frame(xnew) )
    if ( !con )  xnew <- xnew[, -1]
    mu <- xnew %*% be
    est <- mu / sqrt( Rfast::rowsums(mu^2) )
    fit <- NULL
  }

  if ( is.null( colnames(y) ) ) {
    colnames(est) <- colnames(be) <- colnames(seb) <- c("X", "Y", "Z")
  } else  colnames(est) <- colnames(be) <- colnames(seb) <- colnames(y)
  rownames(be) <- rownames(seb) <- c( colnames(x) )

  param <- c(di, rho, psi)
  names(param) = c("Fit measure", "Rho", "Psi")
  theta <- c(the1, the2)
  names(theta) <- c("theta_1", "theta_2")
  list(loglik = -da$value - n * log(2 * pi), param = param, rl = rl,
       theta = theta, beta = be, seb = seb, est = est)
}
