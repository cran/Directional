\name{Density of the spherical Kent and ESAG distributions}
\alias{dkent}
\alias{desag}

\title{
Density of the spherical Kent and ESAG distributions
}

\description{
Density of the spherical Kent and ESAG distributions.
}

\usage{
dkent(y, G, param, logden = FALSE )
desag(y, mu, gam, logden = FALSE)
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
}
\item{mu}{
The mean vector the ESAG distribution, a vector in \eqn{R^3}.
}
\item{gam}{
The two gamma parameters of the ESAG distribution.
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
Gaussian Distribution. Statistics and Computing, 28(3):689--697.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
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
y <- rkent(1000, k = 10, m = m, b = 4)
mod <- kent.mle(y)
dkent( y, G = mod$G, param = mod$param )
}

