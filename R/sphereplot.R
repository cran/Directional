sphereplot <- function(dat, col = NULL, bgcol = "snow") {

  lati <- matrix(seq(90, -90, len = 50)*pi/180, 50, 50, byrow = TRUE)
  longi <- matrix(seq(-180, 180, len = 50)*pi/180, 50, 50)
  x <- cos(lati) * cos(longi)
  y <- cos(lati) * sin(longi)
  z <- sin(lati)

  rgl::open3d()
  ids <- rgl::persp3d(x, y, z, col = "snow", specular = "black", axes = FALSE,
                      box = FALSE, xlab = "", ylab = "", zlab = "",
                      normal_x = x, normal_y = y, normal_z = z, polygon_offset = 1)
  rgl::contourLines3d(ids, list(latitude = function(x, y, z) asin( z / sqrt(x^2 + y^2 + z^2) ) * 180 / pi,
                                longitude = function(x, y, z) atan2(y, x) * 180 / pi) )
  if ( is.null(col) )  col <- rep(2, dim(dat)[1])
  rgl::points3d(dat[, 1], dat[, 2], dat[, 3], col = col, radius = 1, size = 3)

}












