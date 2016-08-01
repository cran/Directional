\name{Random values simulation from a von Mises-Fisher distribution}
\alias{rvmf}
\title{
Random values simulation from a von Mises-Fisher distribution
}

\description{
It generates random vectors following the von Mises-Fisher distribution. The data can be spherical or hyper-spherical.
}

\usage{
rvmf(n, mu, k)
}

\arguments{
\item{n}{
The sample size.
}
\item{mu}{
The mean direction.
}
\item{k}{
The concentration parameter. If k = 0, random values from the spherical uniform will be drwan. Values from a multivariate normal distribution with zero mean vector and the identity matrix as the covariance matrix. Then each vector becomes a unit vector.
}
}

\details{
It uses a rejection smapling as suggested by Andrew Wood (1994). I have mentioned the description of the algorithm as I found it in Dhillon and Sra in 2003.
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

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf}, \link{rfb}, \link{rvonmises}, \link{rmixvmf}
}
}

\examples{
m <- rnorm(4)
m <- m/sqrt(sum(m^2))
x <- rvmf(100, m, 25)
m
vmf(x)
}

\keyword{ Random values simulation }
\keyword{ von Mises-Fisher distribution }
