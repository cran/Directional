###  Simulates from Matrix Fisher distribution on SO(3) using any parameter F (3X3)

rmatrixfisher <- function(n, F) {


   ## Decompose F into UDV' form and find the corresponding between Matrix Fisher parameter F and Bingham parameter A ##
   convertP_SO3_S3<-function(F)  {

    anal <- svd(F)
    U <- anal$u
    V <- anal$v
    D <- diag( anal$d )
    l1 <- 0
    l2 <- 2 * ( D[2] + D[3] )
    l3 <- 2 * ( D[1] + D[3] )
    l4 <- 2 * ( D[1] + D[2] )
    A <- diag( c(l1, l2, l3, l4) )
    list(A = A, U = U, V = V)
   }

   ## convert Bingham samples (x) into Matrix Fisher samples (X) according to John Kent transformation ##
   convert.X <- function(x)  {

    x1 <- x[1]
    x2 <- x[2]
    x3 <- x[3]
    x4 <- x[4]

    X <- matrix( 0, 3, 3 )

    X[1, 1] <- x1^2 + x2^2 - x3^2 - x4^2
    X[2, 1] <- 2 * ( x1 * x4 + x3 * x2 )
    X[3, 1]<-  - 2 * ( x1 * x3 - x2 * x4 )


    X[1, 2] <-  -2 * ( x1 * x4 - x3 * x2 )
    X[2, 2] <- x1^2 + x3^2 - x2^2 - x4^2
    X[3, 2] <- 2 * ( x1 * x2 + x3 * x4 )

    X[1, 3] <- 2 * ( x1 * x3 + x2 * x4 )
    X[2, 3] <-  -2 * ( x1 * x2 - x3 * x4 )
    X[3, 3] <- x1^2 + x4^2 - x2^2 - x3^2

    X
   }


  X <- array( 0, dim = c(3, 3, n) )     ##  will store the Matrix Fisher Samples as an array

    S <- convertP_SO3_S3(F)               ##  convert matrix Fisher parameter F to Bingham parameter A
    A <- S$A
    U <- S$U
    tV <- t( S$V )
    z <- rbingham(n, A)
    for (i in 1:n) {
      X[, , i] <- ( U %*% convert.X( z[i, ] ) ) %*% tV   ## convert Bingham samples to Matrix Fisher samples
    }

    X
}



