################################
#### Cosines to Azimuth/Plunge
#### Eli Amson 10/2017
#### eli.amson1988@gmail.com
#### References: Amson E, Arnold P, Van Heteren AH, Cannoville A, Nyakatura JA. Trabecular architecture in the forelimb epiphyses of extant xenarthrans (Mammalia). Frontiers in Zoology.
################################
cosap <- function(x, y, z) {
  try(if (length(x)!=length(y)|length(x)!=length(z)|length(z)!=length(y)) stop("Dimensions of x, y, and z do not match"))
  A <- NA
  P <- NA
  for (i in 1:length(x)){
    
    ### A[i] (azimuth, theta) ####
    if (z[i] <= 0) {
      if (y[i] < 0) {
        A[i] <- 180 - atan(x[i]/y[i]) * 180 / pi #if z[i] = 0, A[i]  is mod 180 degrees
      } else {
        if ( y[i] > 0 ) {
          A[i] <- 360 * (x[i] > 0) - atan(x[i]/y[i]) * 180 / pi
        } else {
          if (z[i] == 0) {
            A[i] <- 90 + 180 * (x[i] < 0) # same as 270 in this case
          } else {
            A[i] <-  90 + 180 * (x[i] > 0)
          }
        }
      }
    } else {
      if (-y[i] < 0) {
        A[i] <- 180 - atan(x[i]/y[i]) * 180 / pi
      } else {
        if ( -y[i] > 0 ){
          A[i] <- 360 * (-x[i] > 0) - atan(x[i]/y[i]) * 180 / pi
        } else {
          if ( z[i] == 0) {
            A[i] <- 90 + 180 * (x[i] < 0) # same as 270 in this case
          } else {
            A[i] <- 90 + 180 * (x[i] <= 0)
          }
        }
      }
    }
    ##### P[i] (plunge, delta) ####
    P[i] <- asin(-z[i]) * 180 / pi
    if (z[i] > 0)  P[i] <-  - P[i]
  }
  ##### Append results ####
  list(A = A, P = P)
}