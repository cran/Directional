knnreg.tune <- function(y, x, M = 10, A = 10, ncores = 1, res = "eucl",
  type = "euclidean", estim = "arithmetic", mat = NULL, graph = FALSE) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- dim(y)[1]

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
  per <- matrix(nrow = M, ncol = A)

  if (ncores == 1) {

    runtime <- proc.time()

    for (vim in 1:M) {
      ytest <- y[mat[, vim], , drop = FALSE]  ## test set dependent vars
      ytrain <- y[-mat[, vim], , drop = FALSE]  ## train set dependent vars
      xtest <- x[mat[, vim], , drop = FALSE]  ## test set independent vars
      xtrain <- x[-mat[, vim], , drop = FALSE]  ## train set independent vars
      est <- knn.reg(xtest, ytrain, xtrain, k = 1:A, res = res, type = type, estim = estim)
      if (res == "spher") {
        for ( l in 1:A )   per[vim, l] <- 1 - sum( est[[ l ]] * ytest )  / rmat
      } else  for ( l in 1:A )   per[vim, l] <- sum( (est[[ l ]] - ytest)^2 ) / rmat
    }
    mspe <- Rfast::colmeans(per)
    performance <- min(mspe)
    runtime <- proc.time() - runtime

  } else {

    runtime <- proc.time()
    cl <- makePSOCKcluster(ncores)
    registerDoParallel(cl)
    pe <- numeric(A)
    per <- foreach(vim = 1:M, .combine = rbind, .packages = "Rfast",
	         .export = c("knn.reg", "knn", "colOrder", "dista", "colmeans", "colhameans", "rowsums") ) %dopar% {
       ytest <- y[mat[, vim], , drop = FALSE]  ## test set dependent vars
       ytrain <- y[-mat[, vim], drop = FALSE]  ## train set dependent vars
       xtest <- y[mat[, vim], , drop = FALSE]  ## test set independent vars
       xtrain <- x[-mat[, vim], drop = FALSE]  ## train set independent vars
       est <- knn.reg(xtest, ytrain, xtrain, k = 1:A, res = res, type = type, estim = estim)
       if (res == "spher") {
         for ( l in 1:A )   pe[l] <- 1 - sum( est[[ l ]] * ytest )  / rmat
       } else  for ( l in 1:A )   pe[l] <- sum( (est[[ l ]] - ytest)^2 ) / rmat
      return(pe)
    }
    stopCluster(cl)
    runtime <- proc.time() - runtime
    mspe <- Rfast::colmeans(per)
    names(mspe) <- paste("k=", 1:A, sep = "")
    performance <- min(mspe)
    names(performance) <- "mspe"
  }

  if ( graph )  plot(1:A, mspe, xlab = "Nearest neighbours", ylab = "MSPE", type = "b")
  list(crit = mspe, best_k = which.min(mspe) + 1, performance = performance, runtime = runtime)
}




