knnreg.tune <- function(y, x, M = 10, A = 10, ncores = 1, res = "eucl",
  type = "euclidean", estim = "arithmetic", mat = NULL, graph = FALSE) {

  ## y is the multivariate (or univariate) dependent variable
  ## x contains the independent variables(s)
  ## M is the number of folds, set to 10 by default
  ## it is assumed that the training set contains at least 11 observations
  ## A is the highest number of nearest neighbours
  ## res is the type of response, Euclidean ("eucl") or spherical ("spher")
  ## ncores specifies how many cores to use
  ## type is for the distance, Euclidean or Manhattan distance.
  ## which can of course use here.
  ## Type ?dist so see more
  ## estim is either 'arithmetic', 'harmonic'. How to calculate the
  ## estimated value of the Y using the arithmetic mean or the
  ## harmonic mean of the closest observations.
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- dim(y)[1]
  d <- dim(y)[2]

  if ( is.null(mat) ) {
    nu <- sample(1:n, min( n, round(n / M) * M ) )
    ## It may be the case this new nu is not exactly the same
    ## as the one specified by the user
    ## to a matrix a warning message should appear
    options(warn = -1)
    mat <- matrix( nu, ncol = M )
  } else  mat <- mat

  M <- dim(mat)[2]
  rmat <- dim(mat)[1]
  per <- matrix(nrow = M, ncol = A - 1)

  if (type == "spher") {
    runtime <- proc.time()
    ## The k-NN algorithm is calculated R times. For every repetition a
    ## test sample is chosen and its observations are classified
    for (vim in 1:M) {

      ytest <- y[mat[, vim], , drop = FALSE]  ## test set dependent vars
      ytrain <- y[-mat[, vim], drop = FALSE]  ## train set dependent vars
      apo <- tcrossprod(x[mat[, vim], ], x[-mat[, vim], ] )
      apo[ apo >= 1 ] <- 1
      apo <- acos(apo)
	    ina <- as.vector( mat[, -vim] )
      est <- matrix(nrow = rmat, ncol = d)
      bb <- Rfast::colnth(apo, rep(A, rmat) )
      poies <- matrix(0, nrow = A, ncol = rmat)
      disa <- matrix(0, nrow = A, ncol = rmat)
      for ( l in 1:rmat ) {
        sel <- which( apo[, l] <= bb[l] )[1:A]
        disa[, l] <- apo[ sel, l]
        poies[, l] <- ina[ sel ]
      }
      for ( l in 1:c(A - 1) ) {
        k <- l + 1
        if (estim == "arithmetic") {
          for (i in 1:rmat) {
            xa <- cbind(poies[1:k, i], disa[1:k, i])
            qan <- xa[order(xa[, 2]), ]
            a <- qan[1:k, 1]
            yb <- as.matrix( y[a, ] )
            est[i, ] <- Rfast::colmeans( yb )
          }

        } else if (estim == "harmonic") {
          for (i in 1:rmat) {
            xa <- cbind(poies[1:k, i], disa[1:k, i])
            qan <- xa[order(xa[, 2]), ]
            a <- qan[1:k, 1]
            yb <- as.matrix( y[a, ] )
            est[i, ] <- Rfast::colhameans( yb )
          }
        }

        if (res == "spher") {
          est <- est / sqrt( Rfast::rowsums(est^2) )
          per[vim, l] <- 1 - sum( est * ytest )  / rmat
        } else  per[vim, l] <- sum( (est - ytest)^2 ) / rmat
      }
    }
    mspe <- Rfast::colmeans(per)
    bias <- per[ , which.min(mspe)] - Rfast::rowMins(per, value = TRUE) ## apply(per, 1, min)
    estb <- mean( bias )  ## TT estimate of bias
    performance <- c( min(mspe) + estb, estb)
    mspe <- Rfast::colmeans(per)
    runtime <- proc.time() - runtime

  } else {

    if (ncores == 1) {

      runtime <- proc.time()
      for (vim in 1:M) {
        ytest <- y[mat[, vim], , drop = FALSE]  ## test set dependent vars
        ytrain <- y[-mat[, vim], drop = FALSE]  ## train set dependent vars
        apo <- Rfast::dista( x[mat[, vim], , drop = FALSE], x[-mat[, vim], , drop = FALSE], type = type )
	    ina <- as.vector( mat[, -vim] )
        est <- matrix(nrow = rmat, ncol = d)
        bb <- Rfast::colnth(apo, rep(A, rmat) )
        poies <- matrix(0, nrow = A, ncol = rmat)
        disa <- matrix(0, nrow = A, ncol = rmat)
        for ( l in 1:rmat ) {
          sel <- which( apo[, l] <= bb[l] )[1:A]
          disa[, l] <- apo[ sel, l]
          poies[, l] <- ina[ sel ]
        }
        for ( l in 1:c(A - 1) ) {
          k <- l + 1
          if (estim == "arithmetic") {
            for (i in 1:rmat) {
              xa <- cbind(poies[1:k, i], disa[1:k, i])
              qan <- xa[order(xa[, 2]), ]
              a <- qan[1:k, 1]
              yb <- as.matrix( y[a, ] )
              est[i, ] <- Rfast::colmeans( yb )
            }

          } else if (estim == "harmonic") {
            for (i in 1:rmat) {
              xa <- cbind(poies[1:k, i], disa[1:k, i])
              qan <- xa[order(xa[, 2]), ]
              a <- qan[1:k, 1]
              yb <- as.matrix( y[a, ] )
              est[i, ] <- Rfast::colhameans( yb )
            }
          }

          if (res == "spher") {
            est <- est / sqrt( Rfast::rowsums(est^2) )
            per[vim, l] <- 1 - sum( est * ytest )  / rmat
          } else  per[vim, l] <- sum( (est - ytest)^2 ) / rmat
        }
      }
      runtime <- proc.time() - runtime

    } else {

      runtime <- proc.time()
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      pe <- numeric(A - 1)
      per <- foreach(vim = 1:M, .combine = rbind, .packages = "Rfast",
	     .export = c("knn.reg", "dista", "rowsums", "colmeans", "colVars", "colsums") ) %dopar% {
        ytest <- y[mat[, vim], , drop = FALSE]  ## test set dependent vars
        ytrain <- y[-mat[, vim], , drop = FALSE]  ## train set dependent vars
        apo <- Rfast::dista(x[mat[, vim], , drop = FALSE], x[-mat[, vim], , drop = FALSE], type = type )
	    ina <- as.vector( mat[, -vim] )
        est <- matrix(nrow = rmat, ncol = d)
        bb <- Rfast::colnth(apo, rep(A, rmat) )
        poies <- matrix(0, nrow = A, ncol = rmat)
        disa <- matrix(0, nrow = A, ncol = rmat)
        for ( l in 1:rmat ) {
          sel <- which( apo[, l] <= bb[l] )[1:A]
          disa[, l] <- apo[ sel, l]
          poies[, l] <- ina[ sel ]
        }
        for ( l in 1:c(A - 1) ) {
          k <- l + 1
          if (estim == "arithmetic") {
            for (i in 1:rmat) {
              xa <- cbind(poies[1:k, i], disa[1:k, i])
              qan <- xa[order(xa[, 2]), ]
              a <- qan[1:k, 1]
              yb <- as.matrix( y[a, ] )
              est[i, ] <- Rfast::colmeans( yb )
            }

          } else if (estim == "harmonic") {
            for (i in 1:rmat) {
              xa <- cbind(poies[1:k, i], disa[1:k, i])
              qan <- xa[order(xa[, 2]), ]
              a <- qan[1:k, 1]
              yb <- as.matrix( y[a, ] )
              est[i, ] <- Rfast::colhameans( yb )
            }
          }

          if (res == "spher") {
            est <- est / sqrt( Rfast::rowsums(est^2) )
            pe[l] <- 1 - sum( est * ytest )  / rmat
          } else  pe[l] <- sum( (est - ytest)^2 ) / rmat
        }
        return(pe)
      }
      stopCluster(cl)
      runtime <- proc.time() - runtime
    }
    mspe <- Rfast::colmeans(per)
    bias <- per[ ,which.min(mspe)] - Rfast::rowMins(per, value = TRUE)
    estb <- mean( bias )  ## TT estimate of bias
    names(mspe) <- paste("k=", 2:A, sep = "")
    performance <- c( min(mspe) + estb, estb)
    names(performance) <- c("mspe", "estimated bias")
  }

  if ( graph )  plot(2:c(length(mspe) + 1), mspe, xlab = "Nearest neighbours", ylab = "MSPE", type = "b")
  list(crit = mspe, best_k = which.min(mspe) + 1, performance = performance, runtime = runtime)
}




