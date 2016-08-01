# Genearation of three-dimensional rotation matrices using Habeck's algorithm

habeck.rot <- function(F) {

    ana <- svd(F)
    U <- ana$u
    tV <- t(ana$v)
    D <- ana$d
	
    if ( det( U %*% tV ) < 0 ) {
	  U[, 3] <-  - U[, 3]
	  D[3]<-  -D[3]
    }
	
    lamda1 <- D[1]
    lamda2 <- D[2]
    lamda3 <- D[3]

    # Step 2 --- generate Euler angles ###

    # part a - sample alpha and gamma #

    Beta.val <- 0

    kappa_phi <- (lamda1 + lamda2) * { cos(Beta.val/2) } ^ 2
    kappa_shi <- (lamda1 + lamda2) * { sin(Beta.val/2) } ^ 2

    if ( kappa_phi == 0 ) {
      phi <- runif(1, 0, 1)
	  
    } else {
      phi <- rvonmises(1, 0, kappa_phi, rads = TRUE)
    }

    if ( kappa_shi == 0 ) {
      shi <- runif(1, 0, 1)
    
	} else {
      shi <- rvonmises(1, 0, kappa_shi, rads = TRUE) 
    }
  
    u <- rbinom(1, 1, 0.5)
    alpha <- 0.5 * (phi + shi) + pi * u
    gamma <- 0.5 * (phi - shi) + pi * u

    # part b - sample Beta.val 

    kappa_Beta <- (lamda1 + lamda2) * cos(phi) + (lamda1 - lamda2) * cos(shi) + 2 * lamda3
    r <- runif(1, 0, 1)
    x <- 1 + 2 * log( r + (1 - r) * exp(-kappa_Beta) ) / kappa_Beta
    Beta.val <- acos(x)

    # Step 3 --- Return rotation matrix 

    a11 <- cos(alpha) * cos(Beta.val) * cos(gamma) - sin(alpha) * sin(gamma)
    a21 <-  -cos(alpha) * cos(Beta.val) * sin(gamma) - sin(alpha) * cos(gamma)
    a31 <- cos(alpha) * sin(Beta.val)

    a12 <- sin(alpha) * cos(Beta.val) * cos(gamma) + cos(alpha) * sin(gamma)
    a22 <-  -sin(alpha) * cos(Beta.val) * sin(gamma) + cos(alpha) * cos(gamma)
    a32 <- sin(alpha) * sin(Beta.val)

    a13 <-  -sin(Beta.val) * cos(gamma)
    a23 <- sin(Beta.val) * sin(gamma)
    a33 <- cos(Beta.val)

    S <- matrix( c(a11, a21, a31, a12, a22, a32, a13, a23, a33), ncol = 3 )
    U %*% S %*% tV
               
}



