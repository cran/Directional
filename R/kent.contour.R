 ################################
#### Contour plots of the Kent distribution on the sphere
#### Tsagris Michail 06/2014
#### mtsagris@yahoo.gr
################################
kent.contour <- function(k, b) {
  ## k is the concentration parameter
  ## b is the ovalness parameter
  ## b must be less than k/2
  gam <- c(0, k, 0)
  lam <- c(0, -b, b)
  con <- Directional::fb.saddle(gam, lam)[3]
  rho <- sqrt(2)
  x <- seq(-rho, rho, by = 0.01)
  n <- length(x)
  mat1 <- matrix(rep(x^2, n), ncol = n)
  mat2 <- t(mat1)
  z <- sqrt( mat1 + mat2 )
  ind <- ( z^2 < rho^2 )  ## checks if x^2+y^2 < rho^2
  ind[ !ind ] <- NA
  theta <- 2 * asin(0.5 * z)
  xa <- k * cos(theta) + b * (mat1 - mat2) - con
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
                 nlevels = 1000,



                 # Select color function
                  color.palette = colorRampPalette( c( "blue",
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
                 key.axes = {axis(4, col = "black", cex.lab = 1.2)},

                 # Axes labs
                 xlab = "Latitude",
                 ylab = "Longitude",
                 cex.lab = 1.2)
}