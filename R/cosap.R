################################
#### Cosines to Azimuth/Plunge
#### Eli Amson 10/2017
#### eli.amson1988@gmail.com
#### References: Amson E, Arnold P, Van Heteren AH, Cannoville A, Nyakatura JA. Trabecular architecture in the forelimb epiphyses of extant xenarthrans (Mammalia). Frontiers in Zoology.
################################
cosap <- function(x, y, z) {
  if (z <= 0) {
    if (y < 0) {
      A <- 180 - atan(x/y) * 180 / pi  ## if z = 0, A  is mod 180
    } else {
      if ( y > 0 ) {
          A <- 360 * (x > 0) - atan(x/y) * 180 / pi
      } else {
        if (z == 0) {
          A <- 90 + 180 * (x < 0)  ## same as 270 in this case
        } else {
          A <-  90 + 180 * (x > 0)
        }
      }
    }
  } else {
    if (-y < 0) {
      A <- 180 - atan(x/y) * 180 / pi
    } else {
      if ( -y > 0 ){
        A <- 360 * (-x > 0) - atan(x/y) * 180 / pi
      } else {
        if ( z == 0) {
          A <- 90 + 180 * (x < 0)  ## same as 270 in this case
        } else {
          A <- 90 + 180 * (x <= 0)
        }
      }
    }
  }
  ##### P (plunge ,delta)
  P <- asin(-z) * 180 / pi
  if (z > 0)  P <-  - P
  ##### Append results
  list(A = A, P = P)

}


