\name{Density of some (hyper-)spherical distributions}
\alias{dvmf}
\alias{iagd}
\alias{dpurka}
\alias{dspcauchy}
\title{
Density of some (hyper-)spherical distributions
}

\description{
Density of some (hyper-)spherical distributions.
}

\usage{
dvmf(y, mu, k, logden = FALSE )
iagd(y, mu, logden = FALSE)
dpurka(y, theta, a, logden = FALSE)
dspcauchy(y, mu, rho, logden = FALSE)
}

\arguments{
\item{y}{
A matrix or a vector with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{mu}{
The mean direction (unit vector) of the von Mises-Fisher, the IAG, or of the the spherical Cauchy distribution.
}
\item{theta}{
The mean direction (unit vector) of the Purkayastha distribution.
}
\item{k}{
The concentration parameter of the von Mises-Fisher distribution.
}
\item{a}{
The concentration parameter of the Purkayastha distribution.
}
\item{rho}{
The \eqn{\rho} parameter of the spherical Cauchy distribution.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the von Mises-Fisher, of the IAG, of the Purkayastha or of the spherical Cauchy distribution is computed.
}

\value{
A vector with the (log) density values of y.
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through Maximum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70--83

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution. Communications in Statistics-Theory and Methods, 19(6): 1973--1986.

Kato S. and McCullagh P. (2018). Mobius transformation and a Cauchy family on the sphere. arXiv preprint arXiv:1510.07679.

Kato S. & McCullagh P. (2020). Some properties of a Cauchy family on the sphere derived from the Mobius transformations. Bernoulli, 26(4), 3224--3248.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kent.mle}, \link{rkent}, \link{esag.mle}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- rvmf(1000, m = m, k = 10)
dvmf(y, k=10, m)
}

