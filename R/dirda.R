dirda <- function(xnew, x, ina, type = c("vmf", "iag", "esag", "kent", "sc", "pkbd", "purka") ) {

  xnew <- as.matrix(xnew)
  if ( ncol(xnew) == 1 )   xnew <- t(xnew)
  ina <- as.numeric(ina)
  g <- max(ina)
  p <- dim(x)[2]
  mat <- matrix(0, dim(xnew)[1], g)
  est <- matrix(NA, nrow = dim(xnew)[1], 7)

  if ( sum( type == "vmf") == 1 ) {
    mod <- Rfast::multivmf.mle(x, ina, ell = FALSE)
    ki <- mod$ki
    mat <- (p/2 - 1) * log(ki) + ki * tcrossprod(mod$mi, xnew) - log( besselI(ki, p/2 - 1, expon.scaled = TRUE) ) - ki
    est[, 1] <- Rfast::colMaxs(mat)
  }

  if ( sum( type == "iag") == 1 ) {
    for (j in 1:g) {
      mod <- Directional::iag.mle( x[ina == j, ] )
      mat[, j] <- Directional::iagd(xnew, mod$mesi[1, ], logden = TRUE )
      est[, 2] <- Rfast::rowMaxs(mat)
    }
  }

  if ( sum( type == "esag") == 1 ) {
    if ( p == 3 ) {
      for (j in 1:g) {
        mod <- Directional::esag.mle( x[ina == j, ] )
        mat[, j] <- Directional::desag(xnew, mod$mu, mod$gam, logden = TRUE )
      }
    } else {
      mod <- Directional::ESAGd.mle( x[ina == j, ] )
      mat[, j] <- Directional::dESAGd(xnew, mod$mu, mod$gam, logden = TRUE )
    }
    est[, 3] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "kent") == 1 ) {
    if (p == 3) {
        for (j in 1:g) {
          mod <- Directional::kent.mle( x[ina == j, ])
          mat[, j] <- Directional::dkent(xnew, G = mod$G, param = mod$param[1:2], logden = TRUE )
        }
        est[, 4] <- Rfast::rowMaxs(mat)
    }
  }

  if ( sum( type == "sc") == 1 ) {
    for (j in 1:g) {
      mod <- Directional::spcauchy.mle( x[ina == j, ])
      mat[, j] <- Directional::dspcauchy(xnew, mod$mu, mod$rho, logden = TRUE)
    }
    est[, 5] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "sc2") == 1 ) {
    for (j in 1:g) {
      mod <- Directional::spcauchy.mle2( x[ina == j, ])
      mat[, j] <- Directional::dspcauchy(xnew, mod$mu, mod$rho, logden = TRUE)
    }
    est[, 5] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "pkbd") == 1 ) {
    for (j in 1:g) {
      mod <- Directional::pkbd.mle( x[ina == j, ])
      mat[, j] <- Directional::dpkbd(xnew, mod$mu, mod$rho, logden = TRUE)
    }
    est[, 6] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "pkbd2") == 1 ) {
    for (j in 1:g) {
      mod <- Directional::pkbd.mle2( x[ina== j, ])
      mat[, j] <- Directional::dpkbd(xnew, mod$mu, mod$rho, logden = TRUE)
    }
    est[, 6] <- Rfast::rowMaxs(mat)
  }

  if ( sum( type == "purka") == 1 ) {
    for (j in 1:g) {
      mod <- Directional::purka.mle( x[ina == j, ] )
      mat[, j] <- Directional::dpurka(xnew, mod$theta, mod$alpha, logden = TRUE )
    }
    est[, 7] <- Rfast::rowMaxs(mat)
  }

  colnames(est) <- c("vmf", "iag", "esag", "kent", "sc", "pkbd", "purka")
  est
}
