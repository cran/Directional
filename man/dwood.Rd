\name{Density of the Wood bimodal distribution on the sphere}
\alias{dwood}
\title{
Density of the Wood bimodal distribution on the sphere
}

\description{
Density of the Wood bimodal distribution on the sphere.
}

\usage{
dwood(y, param, logden = FALSE)
}

\arguments{
\item{y}{
A matrix containing two columns. The first one is the latitude and the second is the longitude, both expressed in degrees.
}
\item{param}{
A vector with the 5 parameters, in the order they are returned by the \code{\link{wood.mle}} function.
That is, \eqn{(\gamma, \delta, \alpha, \beta, \kappa)}.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the spherical Wood distribution is computed.
}

\value{
A vector with the (log) density values of y.
}

\references{
Wood A.T.A. (1982). A bimodal distribution on the sphere.
Journal of the Royal Statistical Society, Series C, 31(1): 52--58.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dkent}, \link{desag}, \link{wood.mle}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
x <- euclid.inv(x)
mod <- wood.mle(x)
d <- dwood(x, mod$info[, 1])
}

\keyword{ Wood distribution }
\keyword{ Bimodal distribution on the sphere }
