################################
#### Simulating from a Bingham distribution
#### Tsagris Michail 02/2014
#### mtsagris@yahoo.gr
#### References: A new method to simulate the Bingham and related distributions
#### in directional data analysis with applications
#### Kent J.T., Ganeiber A.M. and Mardia K.V. (2013)
#### http://arxiv.org/pdf/1310.8110v1.pdf
#### and
#### Exact Bayesian Inference for the Bingham Distribution
#### C.J. Fallaize and T. Kypraios (2014)
#### http://arxiv.org/pdf/1401.2894v1.pdf
################################

######### Simulation using any symmetric A matrix
rbingham <- function(n, A) {
  Rfast::rbingham(n, A)
}





# rbingham <- function(n, A) {
#   p <- ncol(A)  ## dimensionality of A
#   eig <- eigen(A)
#   lam <- eig$values  ## eigenvalues
#   V <- eig$vectors  ## eigenvectors
#   lam <- lam - lam[p]
#   lam <- lam[-p]
#   ### f.rbing part
#   lam <- sort(lam, decreasing = TRUE)  ## sort the eigenvalues in desceding order
#   nsamp <- 0
#   X <- NULL
#   lam.full <- c(lam, 0)
#   qa <- length(lam.full)
#   mu <- numeric(qa)
#   sigacginv <- 1 + 2 * lam.full
#   SigACG <- sqrt( 1 / ( 1 + 2 * lam.full ) )
#   Ntry <- 0
#
#   while (nsamp < n) {
#     x.samp <- FALSE
#     while ( !x.samp ) {
#       yp <- rnorm(qa, mu, SigACG)
#       y <- yp / sqrt( sum( yp^2 ) )
#       lratio <-  - sum( y^2 * lam.full ) - qa/2 * log(qa) + 0.5 * (qa - 1) + qa/2 * log( sum(y^2 * sigacginv ) )
#       if ( log(runif(1) ) < lratio) {
#         X <- c(X, y)
#         x.samp <- TRUE
#         nsamp <- nsamp + 1
#       }
#       Ntry <- Ntry + 1
#     }
#   }
#
#   x <- matrix(X, byrow = TRUE, ncol = qa)
#   ## the avtry is the estimate of the M in rejection sampling
#   ## 1/M is the probability of acceptance
#   ## the x contains the simulated values
#   tcrossprod(x, V) ## simulated data
# }
