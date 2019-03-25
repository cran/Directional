### MLE of the ESAG
ESAGmle <- function(y, tol = 1e-07) {
  ## y is the spherical data, a matrix with unit vectors
  n <- dim(y)[1]
  I3 <- diag(3)
  z <- t(y)
  nc <- n/2

   mag <- function(param, z, nc){
     m <- param[1:3]
     gam1 <- param[4]
     gam2 <- param[5]
     heta <- sqrt(gam1^2 + gam2^2 + 1) - 1
     m0 <- sqrt( m[2]^2 + m[3]^2 )
     rl <- sum(m^2)
     x1b <- c( -m0^2, m[1] * m[2], m[1] * m[3] ) / m0 / sqrt(rl)
     x2b <- c( 0, -m[3], m[2] )/m0
     T1 <- tcrossprod( x1b )
     T2 <- tcrossprod( x2b )
     T12 <- tcrossprod( x1b, x2b )
     vinv <- I3 + gam1 * ( T1 - T2 ) + gam2 * ( T12 + t(T12) ) + heta * ( T1 + T2 )
     g2 <- colSums( m * z )
     g1 <- colSums( z * crossprod(vinv, z) )
     a <- g2 / sqrt(g1)
     a2 <- a^2
     M2 <- ( 1 + a2 ) * pnorm(a) + a * dnorm(a)
     - 0.5 * sum(a2) + nc * rl + 1.5 * sum( log(g1) ) - sum( log(M2) ) 
   }
  mod <- Rfast::iag.mle(y)
  ini <- as.vector( mod$mesi[1, ] )  ## initial values  
  options(warn = - 1)
  da <- nlm(mag, c( ini, rnorm(2) ), z = z, nc = nc, iterlim = 2000)
  lik1 <-  -da$minimum
  da2 <- nlm(mag, da$estimate, z = z, nc = nc, iterlim = 2000 )
  lik2 <-  -da$minimum
   while ( lik2 - lik1 > tol) {
    lik1 <- lik2
    da <- nlm(mag, da$estimate, z = z, nc = nc, iterlim = 2000 )
    lik2 <-  -da$minimum
   }
  da <- optim(da$estimate, mag, z = z, nc = nc, control = list(maxit = 10000) )
  
  list( mu = da$par[1:3], gam = da$par[4:5], loglik = -da$value - n * log(2 * pi), iag.loglik  = mod$param[2] )
}
