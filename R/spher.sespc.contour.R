spher.sespc.contour <- function(mu, theta, bgcol = "snow", dat = NULL, col = NULL, lat = 50, long = 50) {
  m1 <- mu / sqrt( sum(mu^2) )
  m <- euclid.inv(m1)
  x1 <- seq(m[1] - lat, m[1] + lat, by = 0.1)
  x2 <- seq(m[2] - long, m[2] + long, by = 0.1)
  n <- length(x1)
  wa <- NULL
  for (i in 1:n) {
    w1 <- cbind(x1[i], x2)
    wa <- rbind(wa, w1)
  }
  wa <- Directional::euclid(wa)
  mat <- Directional::dsespc(wa, mu, theta)
  mat <- matrix(mat, nrow = n, byrow = TRUE)
  a <- contourLines(x1, x2, mat, nlevels = 7)

    rgl::open3d()
    lati <- matrix( seq(90, -90, len = 50) * pi / 180, 50, 50, byrow = TRUE )
    longi <- matrix (seq(-180, 180, len = 50) * pi / 180, 50, 50 )
    x <- cos(lati) * cos(longi)
    y <- cos(lati) * sin(longi)
    z <- sin(lati)
    ids <- rgl::persp3d(x, y, z, col = bgcol, axes = FALSE, box = FALSE,
                        xlab = "", ylab = "", zlab = "",
                        normal_x = x, normal_y = y, normal_z = z, polygon_offset = 1)
    rgl::contourLines3d(ids, list(latitude = function(x, y, z) asin( z / sqrt(x^2 + y^2 + z^2) ) * 180 / pi,
                                  longitude = function(x, y, z) atan2(y, x) * 180 / pi) )
    for ( i in 1:length(a) ) {
      y1 <- euclid( cbind( a[[ i ]]$x, a[[ i ]]$y ) )
      rgl::lines3d(y1[, 1], y1[, 2], y1[, 3], col = 4, radius = 1, size = 1)
    }
    rgl::points3d(m1[1], m1[2], m1[3], col = 3, radius = 1, size = 3)
    if ( !is.null(dat) ) {
        if ( is.null(color) )  color <- rep(2, dim(dat)[1])
        rgl::points3d(dat[, 1], dat[, 2], dat[, 3], col = col, radius = 1, size = 3)
    }
}
