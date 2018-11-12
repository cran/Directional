sphereplot <- function(x, col = NULL) {
  
  if ( is.null(col) )  col <- rep(1, dim(x)[1])
  rgl::open3d()
  y1 <- x[, 1]
  y2 <- x[, 2]
  y3 <- x[, 3]
  rgl::points3d(y1, y2, y3, col = col, radius = 1)
  rgl::spheres3d(0, 0, 0, lit = FALSE, color = "white")
  rgl::spheres3d(0, 0, 0, radius = 1, lit = FALSE, color = "black", front = "lines")
}

