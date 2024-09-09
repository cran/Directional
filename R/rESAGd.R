rESAGd <- function(n, mu, gam) {
  d <- length(mu)
  if ( is.null(gam) ) {
    y <- Rfast::matrnorm(n, d) + rep( mu, rep(n, d) )
  } else {	 
    lambda <- .parameter(gam)[[ 1 ]]
    R <- .rotation(gam)
    eigenvector_hat <- .ONB(mu)
    eigenvector <- eigenvector_hat[, -d] %*% R
    P <- cbind(eigenvector, eigenvector_hat[, d])
    V <- P %*% ( t(P) * c(lambda, 1) )
    #V <- .covariance_matrix(mu, lambda, R)
    y <- Rfast::rmvnorm(n, mu, V)
  }
  y / sqrt( Rfast::rowsums(y^2) ) 
}




# Function that map parameter gamma to (lambda, theta, phi)
# @param gamma vector that d = (1+(1+8*(1+g))^0.5)/2 is an integer, where g = length(gamma)
# @return list that contains a vector of lambda, a vector of theta, and a vector of phi
.parameter <- function(gamma) {
  g <- length(gamma)
  d <- 0.5 * ( 1 + sqrt( 1 + 8 * (1 + g) ) )
  if (as.integer(d) != d) { return('incorrect dimensions of the input vector') }
  lambda <- rep(NA, d - 1)
  
  K <- rep(NA, d - 2)  ## creating K for lambda_i=K_i*lambda_{i-1}
  theta <- rep(NA, d - 2)
  dim_phi <- 0.5 * (d - 2) * (d - 3)
  phi <- c()
  
  K[1] <- sqrt( gamma[1]^2 + gamma[2]^2 ) + 1 ## calculate K_1 
  theta[1] <- atan2(gamma[2], gamma[1])  ## calculate theta1
  
  counting <- function(n)  0.5 * n * (n + 1)  ## creating a counting function
  
  if ( d >= 4 ) {
    for ( i in 2:(d - 2) ) {
      c1 <- counting(i)
      c2 <- counting(i + 1) - 1
      group <- gamma[c1:c2] #grouping parameters
      n <- length(group)
      
      K[i] <- sqrt( sum(group^2) ) + 1  ## K
      phi_group <- rep(NA, n - 2)  ## phi
      for ( j in 1:(n - 2) ) {
        if ( sum( group[j:n]^2 ) == 0 ) {
          phi_group[j] <- 0
        } else phi_group[j] <- acos(group[j] / sqrt( sum( group[j:n]^2) ) )  
      }
      phi <- c(phi, phi_group)
      theta[i] <- atan2( group[n], group[n - 1] )
    }  ## end for ( i in 2:(d - 2) ) {
  }  ## end  if ( d >= 4 ) {
  
  K_prod <- prod( K^seq(d - 2 ,1 ) )
  lambda[1] <- (1 / K_prod)^(1 / ( d - 1) )
  for ( j in 2:(d - 1) )  lambda[j] <- K[j - 1] * lambda[j - 1]
 
  list(lambda = lambda, theta = theta, phi = phi)
}


# Function to create Rotation Matrix
# @param gamma a vector that d = (1+(1+8*(1+g))^0.5)/2 is an integer, where g = length(gamma)
# @return Rotation mapping, (d-1)*(d-1) matrix 
.rotation <- function(gamma){
  para <- .parameter(gamma)
  theta <- para[[2]]  
  phi <- para[[3]]  
  d <- length(theta) + 2
  R <- diag(d - 1) #create a (d-1)-by-(d-1) identity matrix
  
  if ( d >= 4 ) {
    for ( m in 1:(d - 3) ) {
      Lo_R <- diag(d - 1)
      Lo_R[1, 1] <- cos( theta[d - m - 1] )
      Lo_R[2, 2] <- Lo_R[1, 1]
      Lo_R[1, 2] <-  -sin( theta[d - m - 1] )
      Lo_R[2, 1] <-  -Lo_R[1, 2]
      
      La_R <- diag(d - 1)
      for ( j in 1:(d - m - 2) ) {
        Rotation <- diag(d - 1) 
        Rotation[j + 1, j + 1] <- cos( phi[1 - j + (d - m - 1) * (d - m - 2 )/2] )
        Rotation[j + 2, j + 2] <- Rotation[j + 1, j + 1]
        Rotation[j + 1, j + 2] <-  -sin( phi[1 - j + (d - m - 1) * (d - m - 2)/2] )
        Rotation[j + 2, j + 1] <-  -Rotation[j + 1, j + 2]
        La_R <- La_R %*% Rotation
      }
      R <- R %*% Lo_R %*% La_R 
    }
  } else  R <- diag(d - 1)
  Lo_R <- diag(d - 1)
  Lo_R[1, 1] <- cos( theta[1] )
  Lo_R[2, 2] <- Lo_R[1, 1]
  Lo_R[1, 2] <-  -sin( theta[1] )
  Lo_R[2, 1] <-  -Lo_R[1, 2]
  R <- R %*% Lo_R
  R
}


# Function to create Variance-Covariance Matrix, V, in ESAG(mu,V)
# @param mu mean direction, non-zero vector with length d
# @lambda eigenvalues of V, vector of positive lambdas with length d-1, generated from function parameter(gamma)
# @R Rotation mapping, (d-1) by (d-1) matrix generated from function rotation(gamma)
# @return d by d matrix
#.covariance_matrix <- function(mu, lambda, R) {
#  d <- length(mu)
#  eigenvector_hat <- .ONB(mu)
#  eigenvector <- eigenvector_hat[, -d] %*% R
#  P <- cbind(eigenvector, eigenvector_hat[, d])
#  V <- P %*% diag( c(lambda, 1) ) %*% t(P)
#  V
#}


# Function to construct orthonormal basis that contains the unit vector mu/|mu|
# @param mu non-zero vector with length d
# @return d*d matrix that columns are the orthonormal basis, and the last column is mu/|mu|
.ONB <- function(mu) {
  d <- length(mu)
  u <- matrix(NA, d, d)
  u[, 1] <- mu 
  u[, 2] <- c( -mu[2], mu[1], rep(0, d - 2) )
  
  for (i in 3:d) {
    u_i <- mu * c( rep(1, i - 1), rep(0, d - i + 1) ) * mu[i]
    u_i[i] <-  -sum( mu[ 1:(i - 1) ]^2 )
    u[, i] <- u_i
  }
  
  #V <- matrix(NA, d, d)
  #std = function(x){ return(x/sqrt(sum(x^2))) }
  #V = apply(u,2,std)
  V <- t( t(u) / sqrt( Rfast::colsums(u^2) ) )
  K <- cbind( V[, 2:d], V[, 1] )
  
  ord <- which( is.na( Rfast::colsums(K) ) ) #to find out columns that are NA (indicating mu_1,...,mu_k = 0 for k>=2) 
  for (i in ord ) { # replace these columns with e_j
    K[, i] <- numeric(d) 
    K[i, i] <- 1
  }
  K
}

