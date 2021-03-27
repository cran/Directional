\name{Simulation of random values from rotationally symmetric distributions}
\alias{rvmf}
\alias{riag}
\title{
Simulation of random values from rotationally symmetric distributions
}

\description{
Simulation of random values from rotationally symmetric distributions. The data can be spherical or hyper-spherical.
}

\usage{
rvmf(n, mu, k)
riag(n, mu)
}

\arguments{
\item{n}{
The sample size.
}
\item{mu}{
A unit vector showing the mean direction for the von Mises-Fisher distribution. The mean vector of the Independent Angular Gaussian distribution. This does not have to be a unit vector.
}
\item{k}{
The concentration parameter of the von Mises-Fisher distribution. If k = 0, random values from the spherical uniform will be drwan.
}
}

\details{
The von Mises-Fisher uses the rejection smapling suggested by Andrew Wood (1994). For the Independent Angular Gaussian, values are generated from a multivariate normal distribution with the given mean vector and the identity matrix as the covariance matrix. Then each vector becomes a unit vector.
}

\value{
A matrix with the simulated data.
}

\references{
Wood A. T. A. (1994). Simulation of the von Mises Fisher distribution. Communications in statistics-simulation and computation, 23(1): 157--164.

Dhillon I. S. & Sra S. (2003). Modeling data using directional distributions. Technical Report TR-03-06, Department of Computer Sciences, The University of Texas at Austin.
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.75.4122&rep=rep1&type=pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.mle}, \link{iag.mle} \link{rfb}, \link{racg}, \link{rvonmises}, \link{rmixvmf}
}
}

\examples{
m <- rnorm(4)
m <- m/sqrt(sum(m^2))
x <- rvmf(100, m, 25)
m
vmf.mle(x)
}

\keyword{ Random values simulation }
\keyword{ von Mises-Fisher distribution }
