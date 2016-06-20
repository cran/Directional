### Code for converting (3 x 3) rotation matrix on SO(3) to unsigned unit quaternion in R4  ###

rot2quat <-function(X)    {

    tr <- X[1, 1] + X[2, 2] + X[3, 3]

    if ( tr > 0 )  {
       s <- sqrt( 1 + tr ) * 2
       x4 <- 0.25 * s
       x1 <- ( X[3, 2] - X[2, 3] ) / s
       x2 <- ( X[1, 3] - X[3, 1] ) / s
       x3 <- ( X[2, 1] - X[1, 2] ) / s
       return( c(x1, x2, x3, x4) )

    }  else  {

       if ( X[1, 1] > X[2, 2]  &  X[1, 1] > X[3, 3] )  {
           s <- sqrt( 1 + X[1, 1] - X[2, 2] - X[3, 3] ) * 2
           x4 <-( X[3, 2] - X[2, 3] ) / s
           x1 <- 0.25 * s
           x2 <- ( X[1, 2] + X[2, 1] ) / s
           x3 <- ( X[1, 3] + X[3, 1] ) / s
           return( c(x1, x2, x3, x4) )

       }  else  {

          if  ( X[2, 2] > X[3, 3] )  {
               s <-sqrt( 1 + X[2, 2] - X[1, 1] - X[3, 3] ) * 2
               x4 <-( X[1, 3] - X[3, 1] ) / s
               x1 <-( X[1, 2] + X[2, 1] ) / s
               x2 <- 0.25 * s
               x3 <- ( X[3, 2] + X[2, 3] ) / s
               return( c(x1, x2, x3, x4) )

          }  else  {

             s <- sqrt( 1 + X[3, 3] - X[1, 1] - X[2, 2] ) * 2
             x4 <- ( X[2, 1] - X[1, 2] ) / s
             x1 <- ( X[1, 3] + X[3, 1] ) / s
             x2 <- ( X[3, 2] + X[2, 3] ) / s
             x3 <- 0.25 * s
             return( c(x1, x2, x3, x4) )
          }

       }


    }

}


