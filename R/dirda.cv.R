dirda.cv <- function(x, ina, folds = NULL, nfolds = 10, k = 2:10, stratified = FALSE,
            seed = FALSE, type = c("vmf", "esag", "sknn", "nsknn"), B = 1000) {

    makefolds <- function(target, nfolds = 10, stratified = TRUE, seed = FALSE) {
    names <- paste("Fold", 1:nfolds)
    runs <- sapply(names, function(x) NULL)
    if (seed)  set.seed(1234)

    if ( !stratified ) {
      options(warn = -1)
      mat <- matrix(sample(length(target)), ncol = nfolds)
      for (i in 1:c(nfolds - 1)) runs[[i]] <- mat[, i]
      names <- prod(dim(mat)) - length(target)
      runs[[nfolds]] <- mat[1:c(nrow(mat) - names), nfolds]
    } else {
      labs <- unique(target)
      run <- list()
      for (i in 1:length(labs)) {
        names <- which(target == labs[i])
        run[[i]] <- sample(names)
      }
      run <- unlist(run)
      for (i in 1:length(target)) {
        k <- i %% nfolds
        if (k == 0)  k <- nfolds
        runs[[k]] <- c(runs[[k]], run[i])
      }
    }
    runs
  }

  if ( is.null(folds) )  folds <- makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)

  est1 <- est2 <- est3 <- est4 <- list()
  for (i in 1:nfolds) {
    est1[[ i ]] <- NA
    est2[[ i ]] <- NA
    est3[[ i ]] <- NA
    est4[[ i ]] <- NA
  }
  per1 <- per2 <- per3 <- per4 <- NULL
  p <- dim(x)[2]

  if ( sum( type == "vmf") == 1 ) {
    per1 <- matrix(0, nfolds, 1)
    colnames(per1) <- "vmf"
    for (i in 1:nfolds) {
      nu <- folds[[ i ]]
      mod <- Rfast::multivmf.mle(x[-nu, ], ina[-nu], ell = FALSE)
      ki <- mod$ki
      mat <- (p/2 - 1) * log(ki) + ki * tcrossprod(mod$mi, x[nu, ]) - log( besselI(ki, p/2 - 1, expon.scaled = TRUE) ) - ki
      est1[[ i ]] <- Rfast::colMaxs(mat)
      per1[i] <- mean(est1[[ i ]] == ina[nu])
    }
  }

  if ( sum( type == "esag") == 1 ) {
    per2 <- matrix(0, nfolds, 1)
    colnames(per2) <- "esag"
    g <- max(ina)
    if (p == 3) {
      for (i in 1:nfolds) {
        nu <- folds[[ i ]]
        mat <- matrix(0, length(nu), g)
        xtrain <- x[-nu, ]
        xtest <- x[nu, ]
        id <- ina[-nu]
        for (j in 1:g) {
          mod <- ESAGmle( xtrain[id == j, ] )
          mat[, j] <- ESAGdensity(xtest, c(mod$mu, mod$gam), logden = TRUE )
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

  if ( sum( type == "sknn" ) == 1 ) {
     per3 <- matrix(0, nfolds, length(k) )
     colnames(per3) <- paste("sknn", k, sep = " ")
     g <- max(ina)
     for (i in 1:nfolds) {
       nu <- folds[[ i ]]
       est3[[ i ]] <- Rfast::dirknn(x[nu, ], x[-nu, ], ina[-nu], k = k, type = "C", parallel = FALSE)
       per3[i, ] <- Rfast::colmeans(est3[[ i ]] == ina[nu])
     }
  }

  if ( sum( type == "nsknn" ) == 1 ) {
     per4 <- matrix(0, nfolds, length(k) )
     colnames(per4) <- paste("nsknn", k, sep = " ")
     g <- max(ina)
     for (i in 1:nfolds) {
       nu <- folds[[ i ]]
       est4[[ i ]] <- Directional::dirknn(x[-nu, ], x[nu, ], k = k, ina[-nu], type = "NS")
       per4[i, ] <- Rfast::colmeans(est4[[ i ]] == ina[nu])
     }
  }

  per <- cbind(per1, per2, per3, per4)
  perf <- Rfast::colmeans(per)
  names(perf) <- colnames(per)
  best <- perf[ which.max(perf) ]
  boot.perf <- NULL
  if (B > 1  & length(perf) > 1 ) {
    est <- cbind(est1[[ 1 ]], est2[[ 1 ]], est3[[ 1 ]], est4[[ 1 ]] )
    for (i in 2:nfolds) {
      est <- rbind(est, cbind(est1[[ i ]], est2[[ i ]], est3[[ i ]], est4[[ i ]] ) )
    }
    diaf <- est - ina[ unlist(folds) ]
    n <- dim(est)[1]
    boot.perf <- numeric(B)
    for (i in 1:B) {
      ind <- sample.int(n, n, TRUE)
      m <- Rfast::colmeans(diaf[ind, ] == 0)
      poio <- which.max(m)
      boot.perf[i] <- mean(diaf[-ind, poio] == 0)
    }
    boot.perf <- mean(boot.perf)
  }
  list(perf = perf, best = best, boot.perf = boot.perf)
}
















