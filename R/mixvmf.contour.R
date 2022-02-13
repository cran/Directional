################################
#### Contour plots of a spherical model based clustering using mixtures
#### of von Mises-Fisher distributions
#### Tsagris Michail 4/2015
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################
mixvmf.contour <- function(u, mod) {
  ## u contains the data in latitude and longitude the first column is
  ## the latitude and the second column is the longitude
  ## mod is a mix.vmf object
  n <- dim(u)[1]  ## sample size
  n1 <- 100
  x1 <- seq(min(u[, 1]) - 5, max(u[, 1]) + 5, length = n1)  ## latitude
  x2 <- seq(min(u[, 2]) - 5, max(u[, 2]) + 5, length = n1)  ## longitude
  mat <- matrix(nrow = n1, ncol = n1)
  mu <- mod$param[, 1:3]  ## mean directions
  tmu <- t(mu)
  k <- mod$param[, 4]  ## concentration parameters
  p <- mod$param[, 5]  ## mixing probabilities
  g <- length(p)  ## how many clusters
  lika <- con <- numeric(g)
  for (l in 1:g)  con[l] <- 0.5 * log(k[l]) - 1.5 * log(2 * pi) - log(besselI(k[l], 0.5, expon.scaled = TRUE)) - k[l]
  for (i in 1:n1) {
    for (j in 1:n1) {
      #x <- c( cos(x1[i]) * cos(x2[j]), cos(x1[i]) * sin(x2[j]), sin(x2[j]) )
      x <- Directional::euclid( c(x1[i], x2[j]) )
      lika <- con + k * ( x %*% tmu )
      can <- sum( p * exp(lika) )
      if (abs(can) < Inf)   mat[i, j] <- can
    }
  }


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
                 plot.axes = { axis(1, col = "black", cex.axis = 1.2);
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