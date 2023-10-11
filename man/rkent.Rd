\name{Simulation of random values from a spherical Kent distribution}
\alias{rkent}
\title{
Simulation of random values from a spherical Kent distribution
}

\description{
Simulation of random values from a spherical Kent distribution.
}

\usage{
rkent(n, k, m, b)
}

\arguments{
  \item{n}{
  The sample size.
}
  \item{k}{
  The concentraion parameter \eqn{\kappa}. It has to be greater than 0.
}
  \item{m}{
  The mean direction (Fisher part).
}
  \item{b}{
  The ovalness parameter, \eqn{\beta}.
}
}

\details{
Random values from a Kent distribution on the sphere are generated. The function generates from a spherical Kent distribution using \code{\link{rfb}} with an arbitrary mean direction and then rotates the data to have the desired mean direction.
}

\value{
A matrix with the simulated data.
}

\references{
Kent J. T., Ganeiber A. M. and Mardia K. V. (2018). A new unified approach for the simulation of a wide class of directional distributions. 
Journal of Computational and Graphical Statistics, 27(2): 291--301.

Kent J.T., Ganeiber A.M. and Mardia K.V. (2013).
A new method to simulate the Bingham and related distributions
in directional data analysis with applications.
http://arxiv.org/pdf/1310.8110v1.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{rfb}, \link{rbingham}, \link{rvmf}, \link{f.rbing}
}
}

\examples{
k <- 15
mu <- rnorm(3)
mu <- mu / sqrt( sum(mu^2) )
A <- diag( c(-5, 0, 5) )
x <- rfb(500, k, mu, A)
kent.mle(x)
y <- rkent(500, k, mu, A[3, 3])
kent.mle(y)
}

\keyword{ Kent distribution }
\keyword{ Simulated data }

