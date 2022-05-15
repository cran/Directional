dirda.cv <- function(x, ina, folds = NULL, nfolds = 10, k = 2:10, stratified = FALSE,
                  type = c("vmf", "iag", "esag", "kent", "knn"),
                  seed = NULL, B = 1000, parallel = FALSE) {

  if ( is.null(folds) )  folds <- Directional::makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
  nfolds <- length(folds)

  est1 <- est2 <- est3 <- est4 <- est5 <- list()
  for (i in 1:nfolds) {
    est1[[ i ]] <- NA
    est2[[ i ]] <- NA
    est3[[ i ]] <- NA
    est4[[ i ]] <- NA
    est5[[ i ]] <- NA
  }
  per1 <- per2 <- per3 <- per4 <- per5 <- NA
  p <- dim(x)[2]

  if ( sum( type == "vmf") == 1 ) {
    per1 <- matrix(0, nfolds, 1)
    colnames(per1) <- "vmf"
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      mod <- Rfast::multivmf.mle(x[-nu, ], ina[-nu], ell = FALSE)
      ki <- mod$ki
      mat <- (p/2 - 1) * log(ki) + ki * tcrossprod(mod$mi, x[nu, , drop = FALSE]) - log( besselI(ki, p/2 - 1, expon.scaled = TRUE) ) - ki
      est1[[ i ]] <- Rfast::colMaxs(mat)
      per1[i] <- mean(est1[[ i ]] == ina[nu])
    }
  }

  if ( sum( type == "iag") == 1 ) {
    per2 <- matrix(0, nfolds, 1)
    colnames(per2) <- "iag"
    g <- max(ina)
    if (p == 3) {
      for (i in 1:nfolds) {
        nu <- folds[[ i ]]
        mat <- matrix(0, length(nu), g)
        xtrain <- x[-nu, ]
        xtest <- x[nu, ]
        id <- ina[-nu]
        for (j in 1:g) {
          mod <- Rfast::iag.mle( xtrain[id == j, ] )
          a <- as.vector(xtest %*% mod$mesi[1, ])
          a2 <- a^2
          pa <- pnorm(a)
          da <- dnorm(a)
          gm <- pa + a2 * pa + a * da
          rl <- mod$param[1]
          mat[, j] <- 0.5 * a2 - 0.5 * rl + log(gm)
        }
        est2[[ i ]] <- Rfast::rowMaxs(mat)
        per2[i] <- mean(est2[[ i ]] == ina[nu])
      }
    } else {
      for (i in 1:nfolds) {
        est2[[ i ]] <- rep(NA, length(folds[[ i ]]) )
        per2[i] <- NA
      }
    }
  }

  if ( sum( type == "esag") == 1 ) {
    per3 <- matrix(0, nfolds, 1)
    colnames(per3) <- "esag"
    g <- max(ina)
    if (p == 3) {
      for (i in 1:nfolds) {
        nu <- folds[[ i ]]
        mat <- matrix(0, length(nu), g)
        xtrain <- x[-nu, ]
        xtest <- x[nu, ]
        id <- ina[-nu]
        for (j in 1:g) {
          mod <- Directional::esag.mle( xtrain[id == j, ] )
          mat[, j] <- Directional::desag(xtest, mod$mu, mod$gam, logden = TRUE )
        }
        est3[[ i ]] <- Rfast::rowMaxs(mat)
        per3[i] <- mean(est3[[ i ]] == ina[nu])
      }
    } else {
      for (i in 1:nfolds) {
        est3[[ i ]] <- rep(NA, length(folds[[ i ]]) )
        per3[i] <- NA
      }
    }
  }

  if ( sum( type == "kent") == 1 ) {

    per4 <- matrix(0, nfolds, 1)
    colnames(per4) <- "kent"
    g <- max(ina)
    if (p == 3) {
      for (i in 1:nfolds) {
        nu <- folds[[ i ]]
        mat <- matrix(0, length(nu), g)
        xtrain <- x[-nu, ]
        xtest <- x[nu, ]
        id <- ina[-nu]
        for (j in 1:g) {
          mod <- Directional::kent.mle( xtrain[id == j, ])
          mat[, j] <- Directional::dkent(xtest, G = mod$G, param = mod$param[1:2] )
        }
        est4[[ i ]] <- Rfast::rowMaxs(mat)
        per4[i] <- mean(est4[[ i ]] == ina[nu])
      }
    } else {
      for (i in 1:nfolds) {
        est4[[ i ]] <- rep(NA, length(folds[[ i ]]) )
        per4[i] <- NA
      }
    }
  }

  if ( sum( type == "knn" ) == 1 ) {
    per5 <- matrix(0, nfolds, length(k) )
    colnames(per5) <- paste("knn", k, sep = " ")
    g <- max(ina)
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      est5[[ i ]] <- Rfast::dirknn(x[nu, , drop = FALSE], x[-nu, ], ina[-nu], k = k, type = "C", parallel = parallel)
      per5[i, ] <- Rfast::colmeans(est5[[ i ]] == ina[nu])
    }
  }

  per <- cbind(per1, per2, per3, per4, per5)
  perf <- Rfast::colmeans(per)
  names(perf) <- colnames(per)
  names(perf)[1:4] <- c("vmf", "iag", "esag", "kent")

  if ( any( !is.na(est5) ) )  {

    est <- est5[[ 1 ]]
    for (i in 2:nfolds)   est <- rbind( est, est5[[ i ]] )
    diaf <- est - ina[unlist(folds)]
    n <- dim(est)[1]
    bf <- numeric(B)

    for (i in 1:B) {
      ind <- sample.int(n, n, TRUE)
      m <- Rfast::colmeans(diaf[ind, ] == 0)
      poio <- which.max(m)
      bf[i] <- mean(diaf[-ind, poio] == 0)
    }
    perf[5] <- mean(bf)
    perf <- perf[1:5]
    names(perf)[5] <- "knn"
  }

  perf

}
