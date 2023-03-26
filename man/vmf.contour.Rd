\name{Contour plots of some rotational symmetric distributions}
\alias{vmf.contour}
\alias{spcauchy.contour}
\alias{purka.contour}
\title{
Contour plots of some rotational symmetric distributions
}

\description{
Contour plots of some rotational symmetric distributions.
}

\usage{
vmf.contour(k)
spcauchy.contour(mu, rho, lat = 50, long = 50)
purka.contour(theta, a, lat = 50, long = 50)
}

\arguments{
\item{k}{
The concentration parameter.
}
\item{mu}{
The mean direction (unit vector) of the von Mises-Fisher, the IAG or the spherical Cauchy distribution.
}
\item{rho}{
The \eqn{\rho} parameter of the spherical Cauchy distribution.
}
\item{theta}{
The median direction for the Purkayastha distribution, a unit vector.
}
\item{a}{
The concentration parameter of the Purkayastha distribution.
}
\item{lat}{
A positive number determing the range of degrees to move left and right from the latitude center. See the example to better understand this argument.
}
\item{long}{
A positive number determing the range of degrees to move up and down from the longitude center. See the example to better understand this argument.
}
}

\details{
The user specifies the concentration parameter only and not the mean direction or data. This is for illustration purposes only. The graph of the von Mises-Fisher distribution will always contain circles, as this distribution is the analogue of a bivariate normal in two dimensions with a zero covariance.
}

\value{
A contour plot of the distribution.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Kato S. and McCullagh P. (2018). Mobius transformation and a Cauchy family on the sphere. arXiv preprint arXiv:1510.07679.

Kato S. & McCullagh P. (2020). Some properties of a Cauchy family on the sphere derived from the Mobius transformations. Bernoulli, 26(4), 3224--3248.

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through Maximum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70--83

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution. Communications in Statistics-Theory and Methods, 19(6): 1973--1986.
}

\author{
Michail Tsagris and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rvmf}, \link{vmf.mle}, \link{vmf.kerncontour}, \link{kent.contour}, \link{sphereplot}
}
}

\examples{
\donttest{
vmf.contour(5)
mu <- colMeans( as.matrix( iris[,1:3] ) )
mu <- mu / sqrt( sum(mu^2) )
spcauchy.contour(mu, 0.7, 30, 30)
spcauchy.contour(mu, 0.7, 60, 60)
}
}

\keyword{ von Mises-Fisher distribution }
\keyword{ Contour plot }
