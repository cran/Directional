knn.reg <- function(xnew, y, x, k = 5, res = "eucl", type = "euclidean", estim = "arithmetic") {
  ## xnew is the new observation
  ## y is the multivariate or univariate dependent variable
  ## x contains the independent variable(s)
  ## k is the number of nearest neighbours to use
  ## res is for the response variable, is it Euclidean or spherical
  ## res = "eucl" for Euclidean, or a real valued one and
  ## res = "spher" for spherical or hyper-spherical
  ## type is for the distance, Euclidean or Manhattan distance.
  ## The function dist()  allows for more distance types
  ## which can of course use here.
  ## Type ?dist so see more
  ## estim is either 'arithmetic', 'harmonic'. How to calculate the
  ## estimated value of the Y using the arithmetic mean or the
  ## harmonic mean of the closest observations

  y <- as.matrix(y)
  d <- ncol(y)  ## dimensions of y
  x <- as.matrix(x)
  p <- ncol(x)  ## dimensions of x
  n <- nrow(y)
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = p)
  nu <- nrow(xnew)

  if ( type == "spher" ) {
    ## calculates distance matrix for (hyper-)spherical data
    x <- x / sqrt( rowSums(x^2) )  ## makes sure x are unit vectors
    xnew <- xnew / sqrt( rowSums(xnew^2) )  ## makes sure x are unit vectors

    dis <- tcrossprod(xnew, x)
    dis[ dis >= 1 ] <- 1
    disa <- acos(dis)

  } else if ( type =="euclidean" || type == "manhattan" ) {
    m <- as.vector( Rfast::colmeans(x) )
    s <- as.vector( Rfast::colVars(x, std = TRUE) )
    x <- ( t(x) - m ) / s ## standardize the independent variables
    x <- t(x)

    if (p == 1) {
      xnew <-as.vector(xnew)
      xnew <- (xnew - m) / s
    } else if ( p > 1 ) {
      xnew <- ( t(xnew) - m ) / s  ## standardize the xnew values
      xnew <- t(xnew)
    }

    ## calculates distance matrix for the Euclidean data
    if (p == 1) {
      x <- as.matrix(x)
      x <- matrix(x, ncol = p)
      xnew <- as.matrix(x)
      xnew <- matrix(xnew, ncol = p)
    }

    disa <- matrix( 0, nu, n )

    if (type == "euclidean") {
       z <- t(x)

       for (i in 1:nu) {
         zz <- z - xnew[i, ]
         disa[i, ] <- sqrt( colSums( zz^2 ) )
       }

    } else if (type == "manhattan") {
      z <- t(x)

      for (i in 1:nu) {
        a <- z - xnew[i, ]
        disa[i, ] <- as.vector( Rfast::colsums( abs(a) ) )
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
      est[i, ] <- as.vector( Rfast::colmeans( yb ) )
    }

  } else if (estim == "harmonic") {
    for (i in 1:nu) {
      xa <- cbind(ina, disa[i, ])
      qan <- xa[order(xa[, 2]), ]
      a <- qan[1:k, 1]
      yb <- as.matrix( y[a, ] )
      est[i, ] <- k / as.vector( Rfast::colsums( yb ) )
    }
  }

  if (res == "spher") {
    y <- y / sqrt( rowSums(y^2) )
  } else  y <- y

  if ( is.null(colnames(y)) ) {
     colnames(est) <- paste("yhat", 1:d, sep = "" )
   } else  colnames(est) <- colnames(y)

  if (d == 1)  est <- as.vector(est)
  est

}
