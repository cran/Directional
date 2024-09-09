pkbd.reg <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-6) {

  x <- model.matrix(y~., as.data.frame(x) )
  if ( !con ) x <- x[, -1, drop = FALSE]
  n <- dim(y)[1]  ;  d <- dim(y)[2] - 1  ;  p <- dim(x)[2]
  H <- matrix(0, (d + 1) * p, (d + 1) * p)
  ind <- matrix(1:( (d + 1) * p), ncol = d + 1 )

  tic <- proc.time()
  be <- Rfast::lmfit(x, y)$be
  mu <- x %*% be
  g2 <- Rfast::rowsums(mu^2)
  a <- Rfast::rowsums(y * mu)
  com <- sqrt(g2 + 1)
  com2 <- com - a
  lik <-  - 0.5 * (d + 1) * sum( log( com2 ) ) - 0.5 * (d - 1) * sum( log(com - 1) - log(g2) )

  muc_y <- mu / com - y
  der <- NULL
  for ( j in 1:(d + 1) ) {
    parta <- Rfast::eachcol.apply(x, muc_y[, j] / com2 )
    partb <- Rfast::eachcol.apply(x, mu[, j] / ( com * (com - 1) ) )
    partc <- Rfast::eachcol.apply(x, mu[, j] / g2)
    der <- c(der, 0.5 * (d + 1) * parta + 0.5 * (d - 1) * partb - (d - 1) * partc)
  }

  for ( i in 1:(d + 1) ) {
    for (j in i:(d + 1) ) {
      if (i == j) {
        ### Jacobian of b1
        a1 <- ( com - mu[, i]^2 / com ) / ( com^2 * com2 )
        upa1 <- crossprod(x, x * a1)
        upa2 <- crossprod(x * muc_y[, i] / com2)

        a2 <- ( com - mu[, i]^2 / com ) / ( com^2 * (com - 1) )
        upb1 <- crossprod(x, x * a2)
        upb2 <- crossprod(x * mu[, i] / ( com *( com - 1) ) )

        upc1 <- crossprod( x / sqrt(g2) )
        upc2 <- crossprod( x * mu[, i] / g2 )

        H[ind[, i], ind[, i]] <- 0.5 * (d + 1) * (upa2 - upa1 ) + 0.5 * (d - 1) * ( upb2 - upb1 ) +
          (d - 1) * (upc1 - 2 * upc2)
      } else {
        ### Jacobian of b12
        a1 <- mu[, i] * mu[, j] / ( com^3 * com2)
        upa1 <- crossprod(x, x * a1)
        upa2 <- crossprod(x * muc_y[, i]/com2, x * muc_y[, j]/com2)

        a2 <- mu[, i] * mu[, j] / ( com^3 * (com - 1) )
        upb1 <- crossprod(x, x * a2)
        upb2 <- crossprod(x * mu[, i] / ( com * (com - 1) ), x * mu[, j] / ( com * (com - 1) ) )

        upc <- crossprod(x * mu[, i] / g2, x * mu[, j] / g2 )

        H[ind[, i], ind[, j]] <- H[ind[, j], ind[, i]] <- 0.5 * (d + 1) * ( upa2 + upa1 ) +
          0.5 * (d - 1) * (upb2 + upb1) - 2 * (d - 1) * upc
      }
    }
  }

  be <- be + solve(H, der)
  mu <- x %*% be
  g2 <- Rfast::rowsums(mu^2)
  a <- Rfast::rowsums(y * mu)
  com <- sqrt(g2 + 1)
  com2 <- com - a
  lik[2] <-  - 0.5 * (d + 1) * sum( log( com2 ) ) - 0.5 * (d - 1) * sum( log(com - 1) - log(g2) )

  k <- 2
  while ( abs(lik[k] - lik[k - 1]) > tol ) {
    k <- k + 1
    muc_y <- mu / com - y
    der <- NULL
    for ( j in 1:(d + 1) ) {
      parta <- Rfast::eachcol.apply(x, muc_y[, j] / com2 )
      partb <- Rfast::eachcol.apply(x, mu[, j] / ( com * (com - 1) ) )
      partc <- Rfast::eachcol.apply(x, mu[, j] / g2)
      der <- c(der, 0.5 * (d + 1) * parta + 0.5 * (d - 1) * partb - (d - 1) * partc)
    }

    for ( i in 1:(d + 1) ) {
      for ( j in i:(d + 1) ) {
        if (i == j) {
          ### Jacobian of b1
          a1 <- ( com - mu[, i]^2 / com ) / ( com^2 * com2 )
          upa1 <- crossprod(x, x * a1)
          upa2 <- crossprod(x * muc_y[, i] / com2)

          a2 <- ( com - mu[, i]^2 / com ) / ( com^2 * (com - 1) )
          upb1 <- crossprod(x, x * a2)
          upb2 <- crossprod(x * mu[, i] / ( com *( com - 1) ) )

          upc1 <- crossprod( x / sqrt(g2) )
          upc2 <- crossprod( x * mu[, i] / g2 )

          H[ind[, i], ind[, i]] <- 0.5 * (d + 1) * (upa2 - upa1 ) + 0.5 * (d - 1) * ( upb2 - upb1 ) +
            (d - 1) * (upc1 - 2 * upc2)
        } else {
          ### Jacobian of b12
          a1 <- mu[, i] * mu[, j] / ( com^3 * com2)
          upa1 <- crossprod(x, x * a1)
          upa2 <- crossprod(x * muc_y[, i]/com2, x * muc_y[, j]/com2)

          a2 <- mu[, i] * mu[, j] / ( com^3 * (com - 1) )
          upb1 <- crossprod(x, x * a2)
          upb2 <- crossprod(x * mu[, i] / ( com * (com - 1) ), x * mu[, j] / ( com * (com - 1) ) )

          upc <- crossprod(x * mu[, i] / g2, x * mu[, j] / g2 )

          H[ind[, i], ind[, j]] <- H[ind[, j], ind[, i]] <- 0.5 * (d + 1) * ( upa2 + upa1 ) +
            0.5 * (d - 1) * (upb2 + upb1) - 2 * (d - 1) * upc
        }  ## end if (i == j)
      }  ##  end for ( j in i:(d + 1) )
    }  ##  for ( i in 1:(d + 1) )

    be <- be + solve(H, der)
    mu <- x %*% be
    g2 <- Rfast::rowsums(mu^2)
    a <- Rfast::rowsums(y * mu)
    com <- sqrt(g2 + 1)
    com2 <- com - a
    lik[k] <-  - 0.5 * (d + 1) * sum( log( com2 ) ) - 0.5 * (d - 1) * sum( log(com - 1) - log(g2) )
  }

  runtime <- proc.time() - tic
  seb <- solve( -H )
  seb <- matrix( sqrt( diag(seb) ), ncol = d + 1)

  if ( is.null(xnew) ) {
    mu <- x %*% be
    g2 <- sqrt( Rfast::rowsums(mu^2) )
    est <- mu / g2
    fit <- sum( y * est )
  } else {
    xnew <- model.matrix(~., data.frame(xnew))
    est <- xnew %*% be
    est <- est / sqrt( Rfast::rowsums(est^2) )
    fit <- NULL
    g2 <- NULL
  }
  loglik <- lik[k] - 0.5 * n * (d - 1) * log(2) +
    n * lgamma(0.5 * (d + 1) ) - n * 0.5 * (d + 1) * log(pi) - n * log(2)

  if ( is.null( colnames(y) ) ) {
    colnames(est) <- colnames(be) <- colnames(seb) <- paste("Y", 1:(d+1), sep = "")
  } else  colnames(est) <- colnames(be) <- colnames(seb) <- colnames(y)
  rownames(be) <- rownames(seb) <- colnames(x)

  list( runtime = runtime, iters = k, loglik = loglik, fit = fit, be = be, seb = seb,
        g2 = g2, est = est )
}
