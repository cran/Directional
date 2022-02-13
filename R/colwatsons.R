colwatsons <- function(u, rads = FALSE) {
  if ( !rads )   u <- u / 180 * pi
  Rfast::colwatsons(u)
}

