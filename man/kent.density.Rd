\name{Density of the spherical Kent and ESAG distributions}
\alias{kent.density}
\alias{ESAGdensity}

\title{
Density of the spherical Kent and ESAG distributions
}

\description{
Density of the spherical Kent and ESAG distributions.
}

\usage{
kent.density(y, G, param, logden = FALSE )
ESAGdensity(y, param, logden = FALSE)
}

\arguments{
\item{y}{
A matrix or a vector with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{G}{
For the Kent distribution only, a 3 x 3 matrix whose first column is the mean direction. The second and third columns are the major and minor axes respectively.
}
\item{param}{
For the Kent distribution a vector with the concentration \eqn{\kappa} and ovalness \eqn{\beta} parameters. The \eqn{\psi} has been absorbed inside the matrix G.

For the ESAG distribution, its parameters, the first three are the mean vector in \eqn{R^3} and the next two are the two gammas.

}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the spherical Kent or spherical ESAG distribution is computed.
}

\value{
A vector with the (log) density values of y.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society,
Series B, 44(1): 71-80.

Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28((3):689--697.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kent.mle}, \link{rkent}, \link{ESAGmle}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- rkent(1000, k = 10, m = m, b = 4)
mod <- kent.mle(y)
kent.density( y, G = mod$G, param = mod$param )
}

