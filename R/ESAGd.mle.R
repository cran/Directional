# Function to find MLEs
# @param Y, dataset with n observations, n*d matrix
# @return list that contains a vector of MLEs, maxmimal log-likelihood value, etc from function optim() with 'BFGS'. It outputs Hessian matrix as well
ESAGd.mle <- function(y, full = FALSE) {

  n <- dim(y)[1]  ;   d <- dim(y)[2]
  p <- d - 1
  start_dim <- 0.5 * (d^2 + d - 2)
  B <- rep(0, start_dim)
  B[1:d] <- 1

  logCd <- log( (2 * pi)^(-0.5 * p) )
  zero <- numeric(d)
  Mp <- matrix(1, n, p)

  joint_log_lik <- function(B, y, d, p, zero, Mp) {
    mu <- B[1:d]
    gamma <- B[-(1:d)]
    R <- .rotation(gamma)
    lambda <- .parameter(gamma)[[ 1 ]]
    eigenvector_hat <- .ONB(mu)
    eigenvector <- eigenvector_hat[, -d] %*% R
    P <- cbind(eigenvector, eigenvector_hat[, d] )
    V <- P %*% ( t(P) * c(lambda, 1) )
    yVy <- Rfast::mahala(y, zero, V)
    a <- as.vector( y %*% mu ) / sqrt(yVy)
    Mp[, 1] <- a * pnorm(a) + dnorm(a)
    Mp[, 2] <- (1 + a^2) * pnorm(a) + a * dnorm(a)
    if ( p >= 3 )  for ( j in 3:p )  Mp[, j] <- a * Mp[, j - 1] + (j - 1) * Mp[, j - 2]
    g2 <- sum(mu^2)
    f <-  - 0.5 * d * log(yVy) + 0.5 * ( a^2 - g2) + log(Mp[, p])
    - sum(f)
  }

  mod <- optim( par = B, joint_log_lik, y = y, d = d, p = p, zero = zero, Mp = Mp,
                method = "BFGS", control = list(maxit = 10000) )
  res <- list( mu = mod$par[1:d], gam = mod$par[-(1:d)], loglik = - mod$value + n * logCd )
  if ( full ) {
    mu <- mod$par[1:d]
    gamma <- mod$par[-(1:d)]
    R <- .rotation(gamma)
    lambda <- .parameter(gamma)[[ 1 ]]
    eigenvector_hat <- .ONB(mu)
    eigenvector <- eigenvector_hat[, -d] %*% R
    P <- cbind(eigenvector, eigenvector_hat[, d] )
    V <- P %*% ( t(P) * c(lambda, 1) )
    res <- list( mu = mod$par[1:d], gam = mod$par[-(1:d)], loglik = - mod$value + n * logCd, vinv = solve(V), lambda = lambda )
  }
  res
}


#
# # Function to find MLEs
# # @param Y, dataset with n observations, n*d matrix
# # @return list that contains a vector of MLEs, maxmimal log-likelihood value, etc from function optim() with 'BFGS'. It outputs Hessian matrix as well
# .ESAG.mle <- function(y){
#
#   start_dim <- ((ncol(y)^2+ncol(y)-2)/2)
#   start <- rep(0,start_dim)
#   start[1:ncol(y)] <- 1
#
#   joint_log_likelihood <- function(B, y = y) {
#     n <- nrow(y)
#     f <- rep(0, n)
#     for (i in 1:n) {
#       f[i] <- .log_likelihood(Bx_i = B, y_i = y[i, ])
#     }
#     sum(f)
#   }
#
#   mod <- optim( par = start, joint_log_likelihood, y= y, method = "BFGS", hessian = TRUE,
#          control = list(fnscale = -1, maxit = 10000) )
# }
#
#
#
# # log likelihood function
# # @param Bx_i, vector with length (d+2)(d-1)/2
# # @param y_i, vector with length d
# # @return log-likelihood, numeric
# .log_likelihood <- function(Bx_i, y_i) {
#   d <- length(y_i)
#   mu <- Bx_i[1:d]
#   gamma <- Bx_i[-(1:d)]
#   R <- .rotation(gamma)
#   lambda <- .parameter(gamma)[[ 1 ]]
#   V_inv <- .covariance_inv_matrix(mu, lambda, R)
#   p <- d-1
#   yVy <- t(y_i) %*% V_inv %*% y_i
#   Cd <- (2 * pi)^(-0.5 * p)
#   a <- as.vector( t(y_i) %*% mu / sqrt(yVy) )
#   Mp <- numeric(p)
#   Mp[1] <- a * pnorm(a) + dnorm(a)
#   Mp[2] <- (1 + a^2) * pnorm(a) + a * dnorm(a)
#   if ( p >= 3 )   for ( j in 3:p )  Mp[j] <- a * Mp[j - 1] + (j - 1) * Mp[j - 2]
#   M_p <- Mp[p]
#   log(Cd) - 0.5 * d * log(yVy) + 0.5 * ( ( t(y_i) %*% mu)^2 / yVy - t(mu) %*% mu ) + log(M_p)
# }
#
#
# # Function to create inverse of Variance-Covariance Matrix, V^-1, in ESAG(mu,V)
# # @param mu, mean direction, non-zero vector with length d
# # @param lambda, eigenvalues of V, vector of positive lambdas with length d-1, generated from function parameter(gamma)
# # @param R, Rotation mapping, (d-1)*(d-1) matrix generated from function rotation(gamma)
# # @return d by d matrix
# .covariance_inv_matrix <- function(mu, lambda, R) {
#   d <- length(mu)
#   eigenvector_hat <- .ONB(mu)
#   eigenvector <- eigenvector_hat[, -d] %*% R
#   P <- cbind(eigenvector, eigenvector_hat[, d] )
#   P %*% diag( c(1/lambda, 1) ) %*% t(P)
# }
