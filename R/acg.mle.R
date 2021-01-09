################################
#### Angular central guassian
#### Tsagris Michail 2/2016
#### mtsagris@yahoo.gr
#### References: Tyler D. E. (1987). Statistical analysis for
#### the angular central Gaussian distribution on the sphere.
#### Biometrika 74(3): 579-589.
################################
acg.mle <- function(x, tol = 1e-07) {
  Rfast::acg.mle(x, tol = tol)
}
