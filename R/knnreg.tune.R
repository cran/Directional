knnreg.tune <- function(y, x, nfolds = 10, A = 10, ncores = 1, res = "eucl",
               estim = "arithmetic", folds = NULL, seed = NULL, graph = FALSE) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- dim(y)[1]

  ina <- 1:n
  if ( is.null(folds) )  folds <- Directional::makefolds(ina, nfolds = nfolds, stratified = FALSE, seed = seed)
  nfolds <- length(folds)

  per <- matrix(nrow = nfolds, ncol = A - 1)

  if (ncores <= 1) {

    runtime <- proc.time()

    for (vim in 1:nfolds) {
      ytest <- y[ folds[[ vim ]], , drop = FALSE]  ## test set dependent vars
      ytrain <- y[ -folds[[ vim ]], , drop = FALSE]  ## train set dependent vars
      xtest <- x[ folds[[ vim ]], , drop = FALSE]  ## test set independent vars
      xtrain <- x[ -folds[[ vim ]], , drop = FALSE]  ## train set independent vars
      est <- Directional::knn.reg(xtest, ytrain, xtrain, k = 2:A, res = res, estim = estim)
      if (res == "spher") {
        for ( l in 1:c(A - 1) )   per[vim, l] <- 1 - mean( est[[ l ]] * ytest )
      } else  for ( l in 1:c(A - 1) )   per[vim, l] <- mean( (est[[ l ]] - ytest)^2 )
    }
    runtime <- proc.time() - runtime
    mspe <- Rfast::colmeans(per)
    names(mspe) <- paste("k=", 2:A, sep = "")
    performance <- min(mspe)

  } else {

    runtime <- proc.time()
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    pe <- numeric(A - 1)
    per <- foreach::foreach(vim = 1:nfolds, .combine = rbind, .packages = "Rfast",
	         .export = c("knn.reg", "knn", "colOrder", "dista", "colmeans", "colhameans", "rowsums") ) %dopar% {
       ytest <- y[ folds[[ vim ]], , drop = FALSE]  ## test set dependent vars
       ytrain <- y[ -folds[[ vim ]], drop = FALSE]  ## train set dependent vars
       xtest <- y[ folds[[ vim ]], , drop = FALSE]  ## test set independent vars
       xtrain <- x[ -folds[[ vim ]], drop = FALSE]  ## train set independent vars
       est <- Directional::knn.reg(xtest, ytrain, xtrain, k = 1:A, res = res, estim = estim)
       if (res == "spher") {
         for ( l in 1:c(A - 1) )   pe[l] <- 1 - mean( est[[ l ]] * ytest )
       } else  for ( l in 1:c(A - 1) )   pe[l] <- mean( (est[[ l ]] - ytest)^2 )
      return(pe)
    }
    parallel::stopCluster(cl)
    runtime <- proc.time() - runtime
    mspe <- Rfast::colmeans(per)
    names(mspe) <- paste("k=", 2:A, sep = "")
    performance <- min(mspe)
  }

  if ( graph ){
    plot(1:A, mspe, xlab = "Nearest neighbours", ylab = "MSPE", type = "b", pch = 16,
         cex.axis = 1.2, cex.lab = 1.2, col = "green", lwd = 2)
    abline(v = 2:A, col = "lightgrey", lty = 2)
    abline(h = seq(min(mspe), max(mspe), length = 10), lty = 2, col = "lightgrey" )
  }
  list(crit = mspe, best_k = which.min(mspe) + 1, performance = performance, runtime = runtime)
}




