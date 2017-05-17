################################
#### Simulating from a von Mises distribution using the algorithm for
#### the von Mises-Fisher distribution
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### References: Wood ATA (1994)
#### Simulation of the von Mises Fisher distribution (Communications in Statistics-Simulation) and
#### Inderjit S. Dhillon and Suvrit Sra (2003)
#### Modeling Data using Directional Distributions (Technical report, The University of Texas at Austin)
################################
rvonmises <- function(n, m, k, rads = TRUE) {
 Rfast::rvonmises(n, m, k, rads)
}

