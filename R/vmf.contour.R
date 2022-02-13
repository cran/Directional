 ################################
#### Contour plots of the von Mises-Fisher distribution on the sphere
#### Tsagris Michail 06/2014
#### mtsagris@yahoo.gr
################################

vmf.contour <- function(k) {
  ## k is the concentration parameter
  rho <- pi/2  ## radius of the circular disc
  x <- seq(-rho, rho, by = 0.01)
  n <- length(x)
  mat <- matrix(rep(x^2, n), ncol = n)
  z <- mat + t(mat)
  theta <- sqrt(z)
  ind <- ( theta < rho )  ## checks if x^2+y^2 < rho^2
  xa <- 0.5 * log(k) + k * cos(theta) - 1.5 * log(2 * pi) - log( besselI(k, 0.5, expon.scaled = TRUE) ) - k
  mat <- exp(xa) * ind

  # Continuous color legend
  # Note that it disappears EVERY BLACK LINE!!!!!!
  # So, for the ones you want, you must do col = "black"
  # For more, see here
  # https://stackoverflow.com/questions/8068366/removing-lines-within-filled-contour-legend
  par(fg = NA)

  # Filled contoure plot in base R
  filled.contour(x, x, mat,

                 # Number of levels
                 # the greater the more interpolate
                 nlevels = 200,



                 # Select color function
                 color.palette =  colorRampPalette( c( "blue",
                                                       "cyan",
                                                       "yellow",
                                                       "red") ),

                 # Adjust axes to points
                 plot.axes = {
                   # # Add points
                   #   points(u[, 1], u[, 2],
                   #          col = "black");

                   # Add contour lines
                   contour(x, x, mat,

                           # Color of contour lines
                           # Otherwise par(fg = NA) will
                           # disappear them...

                           col="black",


                           # Number of levels
                           nlevels = 10,

                           # Size of contour numbers
                           labcex = 0.8,

                           # Width of contour lines
                           lwd = 1.5,

                           add = TRUE) },

                 # Legend tick lines
                 key.axes = {axis(4, col = "black", cex.axis = 1.2)},

                 # Axes labs
                 xlab = "Latitude",
                 ylab = "Longitude",
                 cex.lab = 1.2)
}