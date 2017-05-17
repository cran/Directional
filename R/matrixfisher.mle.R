##      Maximum likelihood estimates for Matrix Fisher parameter F (3X3)
matrixfisher.mle <- function(X) {  
   N <- dim(X)[3]
   Xbar <- apply(X, 1:2, sum)
   svd( Xbar / N ) 
}

