################################
#### Angular central guassian
#### Tsagris Michail 2/2016
#### mtsagris@yahoo.gr
#### References: Tyler D. E. (1987). Statistical analysis for
#### the angular central Gaussian distribution on the sphere.
#### Biometrika 74(3): 579-589.
################################

acg <- function(x) {

  p <- dim(x)[2]
  n <- dim(x)[1]
  mu <- numeric(p)
  lam1 <- cov(x)
  maha <- Rfast::mahala(x, mu, lam1)
  down <- sum( 1 / maha )
  tx <- up <- array( dim = c(n, p, p) )
  for (j in 1:n)   tx[j, , ] <- tcrossprod( x[j, ] )
  up <- tx / maha
  up2 <- colSums(up)
  lam2 <- p * up2 / down
  i <- 2
  while ( sum( abs(lam2 - lam1 ) ) > 1e-10 ) {
    i <- i + 1
    lam1 <- lam2
    maha <- Rfast::mahala(x, mu, lam1)
    down <- sum( 1 / maha )
    up <- tx / maha
    up2 <- colSums(up)
    lam2 <- p * up2 / down
  }

  A <- lam2
  if ( is.null( colnames(x) ) ) {
    colnames(A) <- rownames(A) <- paste("X", 1:p, sep = "")
  } else  colnames(A) <- rownames(A) <- colnames(x)
  list(iter = i, cova = A)
}
