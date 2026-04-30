spher.sespc.contour <- function(mu, theta, bgcol = "snow", dat = NULL, col = NULL, lat = 50, long = 50) {

  m1 <- mu / sqrt( sum(mu^2) )
  m  <- Directional::euclid.inv(m1)

  x1 <- seq(m[1] - lat,  m[1] + lat,  by = 0.1)
  x2 <- seq(m[2] - long, m[2] + long, by = 0.1)

  x1 <- x1[ x1 >= -90  & x1 <= 90  ]
  x2 <- x2[ x2 >= -180 & x2 <= 180 ]

  n1 <- length(x1)
  n2 <- length(x2)

  wa <- NULL
  for (i in 1:n1) {
    w1 <- cbind(x1[i], x2)
    wa <- rbind(wa, w1)
  }
  wa  <- Directional::euclid(wa)
  mat <- Directional::dsespc(wa, mu, theta)
  mat <- matrix(mat, nrow = n1, ncol = n2, byrow = TRUE)

  finite_vals <- mat[ is.finite(mat) ]
  levs <- unique( quantile(finite_vals, probs = seq(0.5, 0.97, length.out = 7)) )
  a <- contourLines(x1, x2, mat, levels = levs)

  rgl::open3d()
  lati  <- matrix( seq(90, -90, len = 50) * pi / 180, 50, 50, byrow = TRUE )
  longi <- matrix( seq(-180, 180, len = 50) * pi / 180, 50, 50 )
  x <- cos(lati) * cos(longi)
  y <- cos(lati) * sin(longi)
  z <- sin(lati)
  ids <- rgl::persp3d(x, y, z, col = bgcol, axes = FALSE, box = FALSE,
                      xlab = "", ylab = "", zlab = "",
                      normal_x = x, normal_y = y, normal_z = z,
                      polygon_offset = 1)
  rgl::contourLines3d(ids,
                      list(latitude  = function(x, y, z) asin(z / sqrt(x^2+y^2+z^2)) * 180/pi,
                           longitude = function(x, y, z) atan2(y, x) * 180/pi))

  for ( i in seq_along(a) ) {
    y1 <- Directional::euclid( cbind(a[[i]]$x, a[[i]]$y) ) * 1.002
    rgl::lines3d(y1[, 1], y1[, 2], y1[, 3], col = 4, lwd = 2)
  }

  rgl::points3d(m1[1], m1[2], m1[3], col = 3, size = 5)
  if ( !is.null(dat) ) {
    if ( is.null(col) )  col <- rep(2, nrow(dat))
    rgl::points3d(dat[, 1], dat[, 2], dat[, 3], col = col, size = 3)
  }
}
