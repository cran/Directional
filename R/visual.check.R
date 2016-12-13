# Visual assessment whether matrix Fisher samples is correctly generated or not
visual.check <- function(x, Fa) {
  n <- dim(x)[3]
  l <- numeric(n)
  for ( i in 1:n )  l[i] <- sum( Fa * x[, , i] )  # Kent Method
  plot( 1:n, l, type = "l", col = "red" )
  abline(h = 20.4, lty = "solid", col = "blue")
  l
}
