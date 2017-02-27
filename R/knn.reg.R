knn.reg <- function(xnew, y, x, k = 5, res = "eucl", type = "euclidean", estim = "arithmetic") {
  ## xnew is the new observation
  ## y is the multivariate or univariate dependent variable
  ## x contains the independent variable(s)
  ## k is the number of nearest neighbours to use
  ## res is for the response variable, is it Euclidean or spherical
  ## res = "eucl" for Euclidean, or a real valued one and
  ## res = "spher" for spherical or hyper-spherical
  ## type is for the distance, Euclidean or Manhattan distance.
  ## which can of course use here.
  ## Type ?dist so see more
  ## estim is either 'arithmetic', 'harmonic'. How to calculate the
  ## estimated value of the Y using the arithmetic mean or the
  ## harmonic mean of the closest observations

  y <- as.matrix(y)
  x <- as.matrix(x)
  d <- dim(y)[2]  ## dimensions of y
  xnew <- as.matrix(xnew)
  p <- dim(x)[2]  ## dimensions of x
  n <- dim(y)[1]
  nu <- dim(xnew)[1]

  if ( type == "spher" ) {
    ## calculates distance matrix for (hyper-)spherical data
    dis <- tcrossprod(xnew, x)
    dis[ dis >= 1 ] <- 1
    disa <- acos(dis)

  } else if ( type =="euclidean" || type == "manhattan" ) {
    m <- Rfast::colmeans(x)
    s <- Rfast::colVars(x, std = TRUE)
    x <- t( ( t(x) - m ) / s )  ## standardize the independent variable
    xnew <- t( ( t(xnew) - m ) / s )  ## standardize the xnew values
    ## calculates distance matrix for the Euclidean data
    disa <- matrix( 0, nu, n )

    if (type == "euclidean") {
       z <- t(x)
       for (i in 1:nu) {
         zz <- z - xnew[i, ]
         disa[i, ] <- sqrt( Rfast::colsums( zz^2 ) )
       }

    } else if (type == "manhattan") {
      z <- t(x)
      for (i in 1:nu) {
        a <- z - xnew[i, ]
        disa[i, ] <- Rfast::colsums( abs(a) )
      }
    }

  }

  ina <- 1:n
  est <- matrix(nrow = nu, ncol = d)

  if (estim == "arithmetic") {
    for (i in 1:nu) {
      xa <- cbind(ina, disa[i, ])
      qan <- xa[order(xa[, 2]), ]
      a <- qan[1:k, 1]
      yb <- as.matrix( y[a, ] )
      est[i, ] <- Rfast::colmeans( yb )
    }

  } else if (estim == "harmonic") {
    for (i in 1:nu) {
      xa <- cbind(ina, disa[i, ])
      qan <- xa[order(xa[, 2]), ]
      a <- qan[1:k, 1]
      yb <- as.matrix( y[a, ] )
      est[i, ] <- k / Rfast::colsums( yb )
    }
  }

  if ( is.null(colnames(y)) ) {
     colnames(est) <- paste("yhat", 1:d, sep = "" )
   } else  colnames(est) <- colnames(y)

  if (d == 1)  est <- as.vector(est)
  est
}
