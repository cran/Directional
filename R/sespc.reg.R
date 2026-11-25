################################################################################
## SESPC regression (Tsagris and Alzeley, 2025, ANZJS 67(1): 77-103)
##
##   Q' y_i ~ SESPC( mu(x_i), theta ),   mu(x_i) = B' x_i,   theta constant.
##
## The SESPC uses the same (mu, theta) parametrisation as the ESAG (theta plays
## the role of gamma), so the rotation Q is needed for the same reason as in
## Paine et al. (2020). The likelihood depends on Q only through the axis
## +-Q e1, because the frame (xi1~, xi2~) is built from e1 x mu: it turns with
## any rotation about e1, and flipping e1 -> -e1 changes only the signs of xi1~
## and xi2~, which cancel in V^{-1}. Hence the search is over axes on a
## hemisphere (Fibonacci grid), not over the whole of O(3).
##
## Every output except theta, rho and psi is given in the ORIGINAL coordinates
## of y. theta is defined relative to the frame induced by the selected Q,
## which is returned so that submodels can be fitted with the same Q.
################################################################################

sespc.reg <- function(y, x, con = TRUE, xnew = NULL, n.axes = 50, Q = NULL) {
  ## y:      n x 3 matrix of unit vectors
  ## x:      covariates (vector, matrix or data.frame)
  ## con:    include an intercept (TRUE) or not
  ## xnew:   new covariate values for prediction (optional)
  ## n.axes: approximate number of axes searched over for Q
  ## Q:      if supplied, use this Q and skip the search (e.g. for submodels)

  X <- model.matrix( ~., data.frame(x) )
  if ( !con )  X <- X[, -1, drop = FALSE]

  if ( is.null(Q) ) {
    wa <- .sespc.hemi.axes(n.axes)
    best <- NULL
    for ( j in 1:nrow(wa) ) {
      Qj <- .sespc.rot.e1(wa[j, ])         ## Qj %*% e1 = wa[j, ]
      mod <- .sespc.reg2(y %*% Qj, X)      ## rows are t(Qj) %*% y_i
      if ( is.null(best) || mod$loglik > best$loglik ) {
        best <- mod
        Q <- Qj
      }
    }
  } else {
    best <- .sespc.reg2(y %*% Q, X)
  }

  .sespc.backtransform(best, Q, y, X, con, xnew)
}


## ---------------------------------------------------------------------------
## Fit for a given (already rotated) response z_i = Q' y_i
## ---------------------------------------------------------------------------
.sespc.reg2 <- function(z, X) {
  n <- dim(z)[1]
  p <- dim(X)[2]
  za <- cbind(z[, 1]^2, z[, 1] * z[, 2], z[, 1] * z[, 3], z[, 1] * z[, 2],
              z[, 2]^2, z[, 2] * z[, 3], z[, 1] * z[, 3], z[, 2] * z[, 3], z[, 3]^2)

  nll <- function(param) {
    be <- matrix( param[-(1:2)], ncol = 3 )
    m <- X %*% be
    g1 <- .sespc.vinv(m, param[1], param[2], za)
    if ( any( !is.finite(g1) ) || any(g1 <= 0) )  return(1e+10)
    a <- Rfast::rowsums(z * m)
    rl <- Rfast::rowsums(m^2)
    E <- g1 * rl + g1 - a^2
    if ( any( !(E > 0) ) )  return(1e+10)
    sqe <- sqrt(E)
    up <- g1 * (rl + 1) * sqe * ( atan2(sqe, -a) - atan2(sqe, a) + pi ) + 2 * a * E
    if ( any( !(up > 0) ) )  return(1e+10)
    - sum( log(up) ) + sum( log(g1) ) + 2 * sum( log(E) )
  }

  ## deterministic starting values: least squares for B, isotropy for theta
  ini <- c( 0, 0, as.vector( .lm.fit(X, z)$coefficients ) )
  val <- optim(ini, nll, method = "BFGS", control = list(maxit = 10000))
  repeat {
    val2 <- optim(val$par, nll, method = "Nelder-Mead", control = list(maxit = 5000))
    val2 <- optim(val2$par, nll, method = "BFGS", control = list(maxit = 10000))
    if ( val$value - val2$value < 1e-06 )  break
    val <- val2
  }
  if ( val2$value < val$value )  val <- val2

  H <- optimHess(val$par, nll)
  covb <- tryCatch( solve(H), error = function(e) matrix(NA, length(val$par), length(val$par)) )
  ## the SESPC normalising constant is (2 pi)^{-2}
  list( loglik = -val$value - 2 * n * log(2 * pi), gam = val$par[1:2],
        be = matrix( val$par[-(1:2)], ncol = 3 ), covb = covb )
}


## ---------------------------------------------------------------------------
## Rotate everything (except theta, rho, psi) back to the original coordinates
## ---------------------------------------------------------------------------
.sespc.backtransform <- function(fit, Q, y, X, con, xnew) {
  p <- dim(X)[2]

  ## mean regression coefficients: mu_orig = Q mu_z  =>  B_orig = B_z Q'
  be <- fit$be %*% t(Q)
  ## covariance: vec(B_orig) = (Q kron I_p) vec(B_z); the thetas are unchanged
  Tm <- diag( 2 + 3 * p )
  Tm[ -(1:2), -(1:2) ] <- kronecker( Q, diag(p) )
  covb <- Tm %*% fit$covb %*% t(Tm)
  se <- suppressWarnings( sqrt( diag(covb) ) )   ## NaN flags a non-PD Hessian
  seb <- matrix( se[-(1:2)], ncol = 3 )
  gam <- fit$gam
  set <- se[1:2]

  ## norms of the fitted mu(x_i) and the measure of fit, in original coordinates
  m <- X %*% be
  rl <- sqrt( Rfast::rowsums(m^2) )
  fitm <- sum(y * m / rl)

  ## rho and psi (relative to the frame induced by Q)
  b <- sqrt( sum(gam^2) )
  rho <- sqrt(b^2 + 1) - b
  psi <- 0.5 * atan2( gam[2], gam[1] )

  nam <- if ( is.null( colnames(y) ) )  c("X", "Y", "Z")  else  colnames(y)
  colnames(be) <- colnames(seb) <- nam
  rownames(be) <- rownames(seb) <- colnames(X)
  names(gam) <- names(set) <- c("Theta_1", "Theta_2")

  ## predicted mean directions for xnew (NULL if xnew is NULL)
  est <- NULL
  if ( !is.null(xnew) ) {
    Xn <- model.matrix( ~., data.frame(xnew) )
    if ( !con )  Xn <- Xn[, -1, drop = FALSE]
    mu <- Xn %*% be
    est <- mu / sqrt( Rfast::rowsums(mu^2) )
    colnames(est) <- nam
  }

  list( loglik = fit$loglik, fit = fitm,
        beta = be, seb = seb, theta = gam, set = set,
        rho = rho, psi = psi, rl = rl, est = est,
        Q = Q, covb = covb )
}


## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
## y' V^{-1} y for each observation (as eq. (15) of Paine et al., 2020)
.sespc.vinv <- function(m, gam1, gam2, za) {
  m0 <- sqrt( m[, 2]^2 + m[, 3]^2 )
  rl <- sqrt( Rfast::rowsums(m^2) )
  x1b <- cbind( -m0^2, m[, 1] * m[, 2], m[, 1] * m[, 3] ) / (m0 * rl)
  x2b <- cbind( 0, -m[, 3], m[, 2] ) / m0
  z12 <- x1b[, 1] * x1b[, 2]  ;  z13 <- x1b[, 1] * x1b[, 3]  ;  z23 <- x1b[, 2] * x1b[, 3]
  tx1 <- cbind(x1b[, 1]^2, z12, z13, z12, x1b[, 2]^2, z23, z13, z23, x1b[, 3]^2)
  w23 <- x2b[, 2] * x2b[, 3]
  tx2 <- cbind(0, 0, 0, 0, x2b[, 2]^2, w23, 0, w23, x2b[, 3]^2)
  vinv <- cbind( x1b[, 1] * x2b, x1b[, 1] * x2b[, 2], 2 * x1b[, 2] * x2b[, 2],
                 x1b[, 3] * x2b[, 2] + x1b[, 2] * x2b[, 3], x1b[, 1] * x2b[, 3],
                 x1b[, 3] * x2b[, 2] + x1b[, 2] * x2b[, 3], 2 * x1b[, 3] * x2b[, 3] )
  heta <- sqrt(gam1^2 + gam2^2 + 1) - 1
  1 + Rfast::rowsums( (gam2 * vinv + gam1 * (tx1 - tx2) + heta * (tx1 + tx2)) * za )
}

## approximately k axes, evenly spread on the hemisphere with first coordinate > 0
.sespc.hemi.axes <- function(k) {
  N <- 2 * k
  i <- seq(0.5, N - 0.5)
  z <- 1 - 2 * i / N
  ph <- pi * (1 + sqrt(5)) * i
  w <- cbind( sqrt(1 - z^2) * cos(ph), sqrt(1 - z^2) * sin(ph), z )
  w[ w[, 1] > 0, , drop = FALSE ]
}

## rotation matrix R with R %*% e1 = q (Rodrigues); stable because q[1] > 0
.sespc.rot.e1 <- function(q) {
  q <- q / sqrt( sum(q^2) )
  v <- c(0, -q[3], q[2])                 ## e1 x q
  cc <- q[1]
  K <- matrix( c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3 )
  diag(3) + K + K %*% K / (1 + cc)
}
