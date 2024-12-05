score.cipc <- function(y, X, rads = TRUE, tol = 1e-6) {

  if ( !rads )  y <- y * pi / 180
  mu <- Directional::cipc.mle(y, rads = rads, tol = tol)$mu
  y <- cbind( cos(y), sin(y) )
  n <- dim(y)[1] 
  H <- matrix(0, 4, 4)
  px <- dim(x)[2]
  stat <- numeric(px)
  
  g2 <- sum(mu^2)
  a <- as.vector(y %*% mu )
  com <- sqrt(g2 + 1)
  com2 <- com - a
  muc_y <- cbind( mu[1] / com - y[, 1], mu[2] / com - y[, 2] ) / com2
  x <- cbind(1, y[, 1] ) ## do note create a x matrix from scratch
  dc1 <- cbind( muc_y[, 1], x[, 2])
  dc2 <- cbind( muc_y[, 2], x[, 2])

  katw <- com^2 * com2
  a1 <- ( com - mu[1]^2 / com ) / katw
  a2 <- ( com - mu[2]^2 / com ) / katw
  a3 <- mu[1] * mu[2] / ( com^3 * com2)

  for ( i in 1:px ) {
  
    x[, 2] <- X[, i]   
    dc1[, 2] <- x[, 2] * muc_y[, 1]
    dc2[, 2] <- x[, 2] * muc_y[, 2]
    der <- c( sum( dc1[, 2] ), sum( dc2[, 2] ) )
    ### Jacobian of b1
    up1 <- crossprod(x, x * a1)
    up2 <- crossprod(dc1)
    H[1:2, 1:2] <- up2 - up1
    ### Jacobian of b2
    up1 <- crossprod(x, x * a2)
    up2 <- crossprod(dc2)
    H[3:4, 3:4] <- up2 - up1
    ### Jacobian of b12
    up1 <- crossprod(x, x * a3)
    up2 <- crossprod(dc1, dc2)
    H[1:2, 3:4] <- H[3:4, 1:2] <- up2 + up1
    
    hinv <- solve(H)[ c(6, 8, 14, 16) ]
    hinv <- matrix(hinv, ncol = 2)
    stat[i] <-  - der %*% hinv %*% der 
  }

  pvalue <- pchisq(stat, 2, lower.tail = FALSE)
  cbind(stat, pvalue)  
}


















