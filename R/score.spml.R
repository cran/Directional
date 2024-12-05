score.spml <- function(y, X, rads = TRUE, tol = 1e-06) {
  ## y is the angular dependent variable
  ## x contains the independent variable(s)
  if ( !rads )  y <- y * pi / 180
  ci <- cos(y)
  si <- sin(y)
  u <- cbind(ci, si)
  n <- dim(u)[1]
  f <-  - 0.5   ;   con <- sqrt(2 * pi) 
  x <- cbind(1, ci) ## do note create a x matrix from scratch
  mu <- as.vector( Rfast::spml.mle(y, tol = tol)$mu )
  tau <- as.vector( u %*% mu )
  ptau <- pnorm(tau)
  rat <- ptau / ( exp(f * tau^2)/con + tau * ptau )
  psit <- tau + rat    
  psit2 <- 2 - tau * rat - rat^2
  com <- Rfast::eachrow( psit * u, mu, oper = "-") 
  
  px <- dim(X)[2]
  stat <- numeric(px)

  for (i in 1:px) {
    x[, 2] <- X[, i]
    der <- crossprod(x, com)
    der <- der[c(2, 4)]
    a11 <- crossprod(x, x * (psit2 * ci^2 - 1) )
    a12 <- crossprod(x, x * (psit2 * ci * si ) )
    a22 <- crossprod(x, x * (psit2 * si^2 - 1 ) )
    der2 <- cbind( rbind(a11, a12), rbind(a12, a22) )
  
    hinv <- solve(der2)[ c(6, 8, 14, 16) ]
    hinv <- matrix(hinv, ncol = 2)
    stat[i] <-  - der %*% hinv %*% der 
  }

  pvalue <- pchisq(stat, 2, lower.tail = FALSE)
  cbind(stat, pvalue) 
}
  
