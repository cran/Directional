knn.reg <- function(xnew, y, x, k = 5, res = "eucl", type = "euclidean", estim = "arithmetic") {
  y <- as.matrix(y)
  x <- as.matrix(x)
  dm <- dim(y)  ## dimensions of y
  d <- dm[2]
  xnew <- as.matrix(xnew)
  nu <- dim(xnew)[1]
  klen <- length(k)
  est <- list()

  if (d == 1  & type == "euclidean") {
    if (estim == "arithmetic") {
      method = "average"
    } else  method = "harmonic"
    g <- Rfast::knn(xnew = xnew, y = y, x = x, k = k, type = "R", method = method)
    for (i in 1:klen)  est[[ i ]] <- g[, i]

  } else {
    if ( type == "spher" ) {
      dis <- tcrossprod(x, xnew)
      dis[ dis >= 1 ] <- 1
      dis[ dis <=  -1 ] <-  -1
      disa <- acos(dis)
    } else   disa <- Rfast::dista(xnew, x, trans = FALSE)

    disa <- Rfast::colOrder(disa)[1:max(k), , drop = FALSE]

	if ( estim == "arithmetic" ) {
      for (j in 1:klen) {
        g <- matrix(nrow = nu, ncol = d)
        knn <- k[j]
        for (i in 1:nu) {
          ind <- disa[1:knn, i]
          g[i, ] <- Rfast::colmeans( y[ind, , drop = FALSE] )
        }
        if (res == "spher") {
          est[[ j ]] <- g / sqrt( Rfast::rowsums(g^2) )
        } else  est[[ j ]] <- g
      }  ## end for (j in klen)

    } else {
      for (j in 1:klen) {
        g <- matrix(nrow = nu, ncol = d)
        knn <- k[j]
        for (i in 1:nu) {
          ind <- disa[1:knn, i]
          g[i, ] <- Rfast::colhameans( y[ind, , drop = FALSE] )
        }
        if (res == "spher") {
          est[[ j ]] <- g / sqrt( Rfast::rowsums(g^2) )
        } else  est[[ j ]] <- g
      }  ## end for (j in klen)
    }   ## end if ( estim == "arithmetic" )
  }  ## end if (d == 1  & type == "euclidean")
  names(est) <- paste("k=", k, sep = "")
  est
}
