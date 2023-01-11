\name{Simulation of random values from a spherical Fisher-Bingham distribution}
\alias{rfb}
\title{
Simulation of random values from a spherical Fisher-Bingham distribution
}

\description{
Simulation of random values from a spherical Fisher-Bingham distribution.
}

\usage{
rfb(n, k, m, A)
}

\arguments{
  \item{n}{
  The sample size.
}
  \item{k}{
  The concentraion parameter (Fisher part). It has to be greater than 0.
}
  \item{m}{
  The mean direction (Fisher part).
}
  \item{A}{
  A symmetric matrix (Bingham part).
}
}

\details{
Random values from a spherical Fisher-Bingham distribution are generated. This functions included the option of simulating from a Kent distribution also.
}

\value{
A matrix with the simulated data.
}

\references{
Kent J.T., Ganeiber A.M. and Mardia K.V. (2013). A new method to simulate the Bingham and related distributions in directional data analysis with applications.
http://arxiv.org/pdf/1310.8110v1.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{rbingham}, \link{rvmf}, \link{rkent}, \link{f.rbing}
}
}

\examples{
k <- 15
mu <- rnorm(3)
mu <- mu / sqrt( sum(mu^2) )
A <- cov(iris[, 1:3])
x <- rfb(50, k, mu, A)
vmf.mle(x) ## fits a von Mises-Fisher distribution to the simulated data
## Next we simulate from a Kent distribution
A <- diag( c(-5, 0, 5) )
n <- 100
x <- rfb(n, k, mu, A) ## data follow a Kent distribution
kent.mle(x) ## fits a Kent distribution
vmf.mle(x) ## fits a von Mises-Fisher distribution
A <- diag( c(5, 0, -5) )
n <- 100
x <- rfb(n, k, mu, A) ## data follow a Kent distribution
kent.mle(x) ## fits a Kent distribution
vmf.mle(x) ## fits a von Mises-Fisher distribution
}

\keyword{ Fisher-Bingham distribution }
\keyword{ Simulated data }
\keyword{ Kent distribution }
