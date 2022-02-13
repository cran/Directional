# Visual assessment whether matrix Fisher samples is correctly generated or not
visual.check <- function(x, Fa) {
  n <- dim(x)[3]
  l <- numeric(n)
  for ( i in 1:n )  l[i] <- sum( Fa * x[, , i] )  # Kent Method
  plot( 1:n, l, type = "l", col = "red", cex.axis = 1.2, cex.lab = 1.2, ylab = "Log prob trace of matrix Fisher dist" )
  abline(h = 20.4, lty = "solid", col = "blue")
  l
}

## slightly faster
#visual.check2 <- function(x, Fa) {
#  n <- dim(x)
#  a <- as.vector( t( apply(x, 2, c) ) )
#  a2 <- a * as.vector( t(Fa) )
#  id <- rep(1:n[3], each = n[1] * n[2])
#  l <- group(a2, id)
#  plot( 1:n[3], l, type = "l", col = "red", cex.axis = 1.2, cex.lab = 1.2, ylab = "Log prob trace of matrix Fisher dist" )
#  abline(h = 20.4, lty = "solid", col = "blue")
#  l
#}
