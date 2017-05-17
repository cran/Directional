## Computing Euler angles from a rotation matrix on SO(3) ##

rot2eul <- function(X) {
                
   v1 <- numeric(3)
   v2 <- numeric(3)              
   x1.13 <- asin( X[3, 1] )
   x2.13 <- pi - x1.13                   
   x1.23 <- atan2( X[3, 2] / cos(x1.13), X[3, 3] / cos(x1.13) )
   x2.23 <- atan2( X[3, 2] / cos(x2.13), X[3, 3] / cos(x2.13) )
   x1.12 <- atan2( X[2, 1] / cos(x1.13), X[1, 1] / cos(x1.13) )
   x2.12 <- atan2( X[2, 1] / cos(x2.13), X[1, 1] / cos(x2.13) )
   v1 <- c( x1.13, x1.23, x1.12 )
   v2 <- c( x2.13, x2.23, x2.12 )                
   list( v1 = v1, v2 = v2 )             
}
        


    


