dirda.cv <- function(x, ina, folds = NULL, nfolds = 10, k = 2:10, stratified = FALSE,
                  type = c("vmf", "iag", "esag", "kent", "sknn", "nsknn"),
                  seed = FALSE, B = 1000, parallel = FALSE) {

  if ( is.null(folds) )  folds <- Directional::makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
  nfolds <- length(folds)

  est1 <- est2 <- est3 <- est4 <- est5 <- est6 <- list()
  for (i in 1:nfolds) {
    est1[[ i ]] <- NA
    est2[[ i ]] <- NA
    est3[[ i ]] <- NA
    est4[[ i ]] <- NA
    est5[[ i ]] <- NA
    est6[[ i ]] <- NA
  }
  per1 <- per2 <- per3 <- per4 <- per5 <- per6 <- NA
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
          mod <- Directional::ESAGmle( xtrain[id == j, ] )
          mat[, j] <- Directional::ESAGdensity(xtest, c(mod$mu, mod$gam), logden = TRUE )
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
          mat[, j] <- Directional::kent.density(xtest, G = mod$G, param = mod$param[1:2] )
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

  if ( sum( type == "sknn" ) == 1 ) {
    per5 <- matrix(0, nfolds, length(k) )
    colnames(per5) <- paste("sknn", k, sep = " ")
    g <- max(ina)
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      est5[[ i ]] <- Rfast::dirknn(x[nu, , drop = FALSE], x[-nu, ], ina[-nu], k = k, type = "C", parallel = parallel)
      per5[i, ] <- Rfast::colmeans(est5[[ i ]] == ina[nu])
    }
  }

  if ( sum( type == "nsknn" ) == 1 ) {
    per6 <- matrix(0, nfolds, length(k) )
    colnames(per6) <- paste("nsknn", k, sep = " ")
    g <- max(ina)
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      est6[[ i ]] <- Directional::dirknn(x[-nu, ], x[nu, , drop = FALSE], k = k, ina[-nu], type = "NS")
      per6[i, ] <- Rfast::colmeans(est6[[ i ]] == ina[nu])
    }
  }

  per <- cbind(per1, per2, per3, per4, per5, per6)
  perf <- Rfast::colmeans(per)
  names(perf) <- colnames(per)
  names(perf)[1:4] <- c("vmf", "iag", "esag", "kent")

  if ( any( !is.na(est5) ) &  any( !is.na(est6) ) ) {

    est1 <- est5[[ 1 ]]
    est2 <- est6[[ 1 ]]
    for (i in 2:nfolds) {
      est1 <- rbind( est1, est5[[ i ]] )
      est2 <- rbind( est2, est6[[ i ]] )
    }
    diaf1 <- est1 - ina[unlist(folds)]
    diaf2 <- est2 - ina[unlist(folds)]
    n <- dim(est1)[1]
    bf1 <- bf2 <- numeric(B)

    for (i in 1:B) {
      ind <- sample.int(n, n, TRUE)
      m <- Rfast::colmeans(diaf1[ind, ] == 0)
      poio <- which.max(m)
      bf1[i] <- mean(diaf1[-ind, poio] == 0)
      m <- Rfast::colmeans(diaf2[ind, ] == 0)
      poio <- which.max(m)
      bf2[i] <- mean(diaf2[-ind, poio] == 0)
    }
    perf[5] <- mean(bf1)
    perf[6] <- mean(bf2)
    perf <- perf[1:6]
    names(perf)[5:6] <- c("sknn", "nsknn")

  } else if ( any( !is.na(est5) ) )  {

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
    perf[6] <- NA
    perf <- perf[1:6]
    names(perf)[5] <- c("sknn")
    names(perf)[6] <- c("nsknn")

  } else if ( any( !is.na(est6) ) ){

    est <- est6[[ 1 ]]
    for (i in 2:nfolds)   est <- rbind( est, est6[[ i ]] )
    diaf <- est - ina[unlist(folds)]
    n <- dim(est)[1]
    bf <- numeric(B)

    for (i in 1:B) {
      ind <- sample.int(n, n, TRUE)
      m <- Rfast::colmeans(diaf[ind, ] == 0)
      poio <- which.max(m)
      bf[i] <- mean(diaf[-ind, poio] == 0)
    }
    perf[6] <- mean(bf)
    perf[5] <- NA
    perf <- perf[1:6]
    names(perf)[5] <- c("sknn")
    names(perf)[6] <- c("nsknn")

  }

  perf

}
