spcauchy.contour <- function(mu, rho, lat = 50, long = 50) {
  m <- euclid.inv(mu)
  x1 <- x2 <- seq(0, 180, by = 2)
  x1 <- seq(m[1] - lat, m[1] + lat)
  x2 <- seq(m[2] - long, m[2] + long)
  n <- length(x1)

  wa <- NULL
  for ( i in 1:n ) {
    w1 <- cbind(x1[i], x2)
    wa <- rbind(wa, w1)
  }

  wa <- Directional::euclid(wa)
  mat <- Directional::dspcauchy(wa, mu, rho)
  mat <- matrix(mat, nrow = n, byrow = TRUE)

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
                 # Axes labs
                 xlab = "Latitude",
                 ylab = "Longitude",
                 cex.lab = 1.2)
}
