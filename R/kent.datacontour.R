 ################################
#### Contour plots of the Kent distribution on the sphere
#### with data appearing
#### Tsagris Michail 06/2014
#### mtsagris@yahoo.gr
################################
kent.datacontour <- function(x) {
  ## u contains the data in latitude and longitude
  ## the first column is the latitude and the
  ## second column is the longitude

  ## if u are eucliean coordinates turn them into
  ## latitude and longitude
  dm <- dim(x)
  if ( dm[2] == 3 ) {
    u <- Directional::euclid.inv(x)
  } else if ( dm[2] == 2 ) {
    u <- x
    x <- Directional::euclid(x) ## Euclidean coordinates used by Kent (1982)
  }

  n <- dm[1]  ## sample size
  a <- Directional::kent.mle(x) ## MLE estimation of the Kent distribution
  G <- a$G ## G matrix, the mean direction and the major-minor axes
  k <- a$param[1] ## kappa, concentration parameter
  b <- a$param[2] ## beta, ovalness
  gam <- c(0, k, 0)
  lam <- c(0, -b, b)
  ckb <- Directional::fb.saddle(gam, lam)[3] ## logarithm of the normalising constant
  n <- 100
  x1 <- seq(min(u[, 1]) - 5, max(u[, 1]) + 5, length = n)  ## latitude
  x2 <- seq(min(u[, 2]) - 5, max(u[, 2]) + 5, length = n)  ## longitude

  wa <- NULL
  for ( i in 1:n )  wa <- rbind(wa, Directional::euclid( cbind(x1[i], x2) ) )
  can <- Directional::dkent(wa, a$G, a$param, logden = FALSE)
  mat <- matrix(can, byrow = TRUE, nrow = n, ncol = n)

  # Continuous color legend
  # Note that it disappears EVERY BLACK LINE!!!!!!
  # So, for the ones you want, you must do col = "black"
  # For more, see here
  # https://stackoverflow.com/questions/8068366/removing-lines-within-filled-contour-legend

  par(fg = NA)

  # Filled contoure plot in base R
  filled.contour(x1, x2, mat,

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
                    # Plot axes
                      axis(1, col = "black", cex.axis = 1.2);
                      axis(2, col = "black", cex.axis = 1.2);

                    # Add points
                      points(u[, 1], u[, 2],
                             col = "black");

                   # Add contour lines
                   contour(x1, x2, mat,

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

                 cex.lab = 1.2,

                 # Axes labs
                 xlab = "Latitude",
                 ylab = "Longitude")


}
