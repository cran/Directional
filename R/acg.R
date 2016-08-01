################################
#### Angular central guassian
#### Tsagris Michail 2/2016
#### mtsagris@yahoo.gr
#### References: Tyler D. E. (1987). Statistical analysis for
#### the angular central Gaussian distribution on the sphere.
#### Biometrika 74(3): 579-589.
################################

acg <- function(x) {

  x <- as.matrix(x)
  x <- x /sqrt( rowSums(x^2) )
  p <- ncol(x)
  n <- nrow(x)
  mu <- numeric(p)

  lam1 <- cov(x)
  maha <- Rfast::mahala(x, mu, lam1)
  down <- sum( 1 / maha )
  tx <- up <- array( dim = c(n, p, p) )

  for (j in 1:n) {
    tx[j, , ] <- crossprod( t( x[j, ] ) )
  }
  up <- tx / maha


  up2 <- apply(up, 2:3, sum)
  lam2 <- p * up2 / down
  i <- 2

  y <- t(x)

  while ( sum( abs(lam2 - lam1 ) ) > 1e-10 ) {
    i <- i + 1
    lam1 <- lam2
    sa <- solve(lam1)
    maha <- as.vector( Rfast::colsums( y * crossprod(sa, y) ) )
    down <- sum( 1 / maha )

    up <- tx / maha

    up2 <- apply(up, 2:3, sum)
    lam2 <- p * up2 / down
  }

  A <- lam2

  if ( is.null( colnames(x) ) ) {
    colnames(A) <- rownames(A) <- paste("X", 1:p, sep = "")
  } else  colnames(A) <- rownames(A) <- colnames(x)
  list(iter = i, cova = A)

}
