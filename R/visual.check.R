# Visual assessment whether matrix Fisher samples is correctly generated or not

visual.check <- function(n, F) {

    l1 <- numeric(n)
    l2 <- numeric(n)

    for ( i in 1:n ) {

	  X1 <- habeck.rot(F)                   # Habeck Method
      X2 <- rmatrixfisher(1, F)              # Kent Method
      l1[i] <- sum( diag( t(F) %*% X1 ) )
      l2[i] <- sum( diag( t(F) %*% X2[, , 1] ) )

    }

    list(l1 = l1, l2 = l2)

}



