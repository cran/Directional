esag.nr <- function(y, tol = 1e-6, maxiters = 100) {
  n <- dim(y)[1]
  m <- Rfast::colmeans(y)
  Rbar <- sqrt( sum(m^2) )
  ## rotate so that the sample mean direction is e3; this keeps the frame (14)
  ## away from its singularity at mu2 = mu3 = 0
  e3 <- c(0, 0, 1)  ;  md <- m / Rbar
  v <- md - e3
  rot <- if ( sum(v^2) < 1e-12 )  diag(3)  else  diag(3) - 2 * tcrossprod(v) / sum(v^2)
  rot <- diag( c(1, -1, 1) ) %*% rot      ## Householder reflection -> proper rotation
  yr <- y %*% t(rot)                      ## rot %*% md = e3
  ## IAG: vMF-type start for ||mu||, then Newton on mu with gam = 0
  mod <- Rfast::iag.mle(yr)
  mu <- mod$mesi[1, ]
  ## ESAG start for gam from the tangent-plane covariance (Proposition 2):
  ## cov(u) is proportional to (I + A)^{-1}, A as in (18), with det(I + A) = 1
  u <- yr %*% .esag.frame(mu)
  P <- solve( crossprod(u) )
  P <- P / sqrt( det(P) )
  gam <- c( (P[1, 1] - P[2, 2]) / 2, P[1, 2] )
  ## and, from the same approximation, ||mu||^2 * u'(I + A)u ~ chi^2_2 with P = I + A
  mu <- mu / sqrt( sum(mu^2) ) * sqrt( 2 / mean( Rfast::rowsums( (u %*% P) * u ) ) )
  der <- .esag.derivs(yr, mu, gam)
  th <- c(mu, gam)
  j <- 0  ;  dif <- 1
  while ( dif > tol  &  j < maxiters ) {
    j <- j + 1
    th <- th + .nr.step(der$hess, der$score)
    lik0 <- der$loglik
    der <- .esag.derivs(yr, th[1:3], th[4:5])
    dif <- abs(der$loglik - lik0)
  }
  ## back to the original coordinates: mu rotates, gam is re-read from V^{-1}
  ## through the frame of the rotated-back mu (gam1 = (a11 - a22)/2, gam2 = a12)
  mu <- as.vector( crossprod(rot, th[1:3]) )
  X <- .esag.frame(mu)
  B <- t(X) %*% ( t(rot) %*% .esag.vinv(th[1:3], th[4:5]) %*% rot ) %*% X
  gam <- c( (B[1, 1] - B[2, 2]) / 2, B[1, 2] )
  nam <- c("mu1", "mu2", "mu3", "gam1", "gam2")
  names(mu) <- nam[1:3]  ;  names(gam) <- nam[4:5]
  list( mu = mu, gam = gam, loglik = der$loglik, iters = j, iag.loglik = mod$param[2] )
}


## psi1 = M1/phi, psi2 = M2/phi and R = Phi/phi, computed stably
.esag.psi <- function(a) {
  R <- exp( pnorm(a, log.p = TRUE) - dnorm(a, log = TRUE) )
  p1 <- a * R + 1
  p2 <- (1 + a^2) * R + a
  neg <- a < -6
  if ( any(neg) ) {  ## asymptotic series avoids the cancellation for alpha << 0
    b <- -a[neg]
    k <- 0:12
    c1 <- (-0.5)^k * factorial(1 + 2 * k) / factorial(k)
    c2 <- (-0.5)^k * factorial(2 + 2 * k) / factorial(k)
    p1[neg] <- as.vector( outer(b, -(2 + 2 * k), "^") %*% c1 )
    p2[neg] <- as.vector( outer(b, -(3 + 2 * k), "^") %*% c2 )
  }
  list(R = R, p1 = p1, p2 = p2)
}

## log-likelihood, score and observed Hessian at (mu, gam)
.esag.derivs <- function(y, mu, gam, hessian = TRUE) {
  n <- dim(y)[1]
  y1 <- y[, 1]  ;  y2 <- y[, 2]  ;  y3 <- y[, 3]
  m2 <- sum(mu^2)  ;  m <- sqrt(m2)
  mu0 <- sqrt( mu[2]^2 + mu[3]^2 )
  ## frame xi1, xi2 as in (14):  u1 = b / D, u2 = a / mu0
  a <- mu[2] * y3 - mu[3] * y2
  cc <- mu[2] * y2 + mu[3] * y3
  b <- -mu0^2 * y1 + mu[1] * cc
  D <- mu0 * m
  u1 <- b / D  ;  u2 <- a / mu0
  ## gradients of u1, u2 w.r.t. mu (n x 3)
  g0 <- c(0, mu[2], mu[3]) / mu0                  ## d mu0 / d mu
  gD <- g0 * m + mu0 * mu / m                     ## d D / d mu
  Ga <- cbind(0, y3, -y2)                         ## d a / d mu
  Gb <- cbind(cc, -2 * mu[2] * y1 + mu[1] * y2, -2 * mu[3] * y1 + mu[1] * y3)
  G1 <- Gb / D - outer(b, gD) / D^2
  G2 <- Ga / mu0 - outer(a, g0) / mu0^2
  ## q = 1 + u'Au
  r <- sqrt( 1 + sum(gam^2) )
  A11 <- gam[1] + r - 1  ;  A22 <- -gam[1] + r - 1  ;  A12 <- gam[2]
  v1 <- A11 * u1 + A12 * u2  ;  v2 <- A12 * u1 + A22 * u2
  q <- 1 + u1 * v1 + u2 * v2
  s2 <- u1^2 + u2^2
  qmu <- 2 * (G1 * v1 + G2 * v2)
  qgam <- cbind( u1^2 - u2^2 + gam[1] / r * s2, 2 * u1 * u2 + gam[2] / r * s2 )
  ## the scalar kernel h(t, q)
  tt <- as.vector(y %*% mu)
  sq <- sqrt(q)
  al <- tt / sq
  ps <- .esag.psi(al)
  s <- al + 2 * ps$p1 / ps$p2            ## s(alpha) = alpha + 2 M1/M2 = d log psi2 / d alpha
  ds <- 2 * ps$R / ps$p2 - 4 * (ps$p1 / ps$p2)^2 + 1   ## s'(alpha)
  ht <- s / sq
  hq <- -(3 + al * s) / (2 * q)
  loglik <- sum( -1.5 * log(q) + log(ps$p2) ) - 0.5 * n * m2 - 1.5 * n * log(2 * pi)
  smu <- as.vector( crossprod(y, ht) + crossprod(qmu, hq) ) - n * mu
  sgam <- as.vector( crossprod(qgam, hq) )
  score <- c(smu, sgam)
  if ( !hessian )  return( list(loglik = loglik, score = score) )

  htt <- ds / q
  htq <- -(s + al * ds) / (2 * q * sq)
  hqq <- (6 + 3 * al * s + al^2 * ds) / (4 * q^2)
  ## sum_i w_i * Hess(u1_i) and sum_i w_i * Hess(u2_i), using linearity in y
  H0 <- ( diag( c(0, 1, 1) ) - tcrossprod(g0) ) / mu0
  HD <- H0 * m + ( tcrossprod(g0, mu) + tcrossprod(mu, g0) ) / m +
    mu0 * ( diag(3) / m - tcrossprod(mu) / m^3 )
  sumH1 <- function(w) {
    sy <- as.vector( crossprod(y, w) )
    Hb <- matrix( c(0, sy[2], sy[3], sy[2], -2 * sy[1], 0, sy[3], 0, -2 * sy[1]), 3, 3 )
    Sb <- as.vector( crossprod(Gb, w) )
    Hb / D - ( tcrossprod(Sb, gD) + tcrossprod(gD, Sb) ) / D^2 +
      sum(w * b) * ( 2 * tcrossprod(gD) / D^3 - HD / D^2 )
  }
  sumH2 <- function(w) {
    Sa <- as.vector( crossprod(Ga, w) )
    -( tcrossprod(Sa, g0) + tcrossprod(g0, Sa) ) / mu0^2 +
      sum(w * a) * ( 2 * tcrossprod(g0) / mu0^3 - H0 / mu0^2 )
  }
  ## sum_i hq_i * d2 q_i / d mu d mu'
  Qmm <- 2 * ( crossprod(G1, G1 * (A11 * hq)) + crossprod(G1, G2 * (A12 * hq)) +
                 crossprod(G2, G1 * (A12 * hq)) + crossprod(G2, G2 * (A22 * hq)) +
                 sumH1(hq * v1) + sumH2(hq * v2) )
  Ymq <- crossprod(y, qmu * htq)
  Hmm <- crossprod(y, y * htt) + Ymq + t(Ymq) + crossprod(qmu, qmu * hqq) + Qmm - n * diag(3)
  ## cross block
  w1 <- (1 + gam[1] / r) * u1  ;  w2 <- (-1 + gam[1] / r) * u2       ## A_gam1 u
  z1 <- gam[2] / r * u1 + u2   ;  z2 <- u1 + gam[2] / r * u2         ## A_gam2 u
  qmg <- cbind( 2 * crossprod(G1 * w1 + G2 * w2, hq), 2 * crossprod(G1 * z1 + G2 * z2, hq) )
  Hmg <- crossprod(y, qgam * htq) + crossprod(qmu, qgam * hqq) + qmg
  ## gamma block
  Hgg <- crossprod(qgam, qgam * hqq) + sum(hq * s2) * ( diag(2) / r - tcrossprod(gam) / r^3 )
  hess <- rbind( cbind(Hmm, Hmg), cbind(t(Hmg), Hgg) )
  list(loglik = loglik, score = score, hess = hess)
}

## the frame (14) and V^{-1} (18)
.esag.frame <- function(mu) {
  m <- sqrt( sum(mu^2) )  ;  mu0 <- sqrt( mu[2]^2 + mu[3]^2 )
  cbind( c(-mu0^2, mu[1] * mu[2], mu[1] * mu[3]) / (mu0 * m), c(0, -mu[3], mu[2]) / mu0 )
}
.esag.vinv <- function(mu, gam) {
  X <- .esag.frame(mu)
  r <- sqrt( 1 + sum(gam^2) )
  A <- matrix( c(gam[1] + r - 1, gam[2], gam[2], -gam[1] + r - 1), 2, 2 )
  diag(3) + X %*% A %*% t(X)
}

## Newton step; if the Hessian is not negative definite its eigenvalues are
## replaced by -|eigenvalues| (modified Newton, still no line search)
.nr.step <- function(H, s) {
  e <- eigen(H, symmetric = TRUE)
  if ( e$values[1] >= 0 ) {
    return( as.vector( e$vectors %*% ( crossprod(e$vectors, s) / abs(e$values) ) ) )
  }
  -solve(H, s)
}



