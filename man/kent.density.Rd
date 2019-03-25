\name{Density of the von Mises-Fisher and the spherical Kent distributions}
\alias{kent.density}
\alias{vmf.density}
\title{
Density of the spherical Kent distribution
}

\description{
Density of the spherical Kent distribution.
}

\usage{
kent.density(y, G, param, logden = FALSE )
vmf.density(y, k, mu, logden = FALSE)
}

\arguments{
\item{y}{
A matrix or a vector with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{G}{
A 3 x 3 matrix whose first column is the mean direction. The second and third columns are the major and minor axes respectively.
}
\item{param}{
A vector with the concentration \eqn{\kappa} and ovalness \eqn{\beta} parameters. The \eqn{\psi} has been absorbed inside the matrix G.
}
\item{k}{
The concentration parameter of the von Mises-Fisher distribution.
}
]\item{mu}{
The mean direction of the von Mises-Fisher distribution.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the von Mises-Fisher or of the spherical Kent distribution is computed.
}

\value{
A vector with the density values of y.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society,
Series B, 44(1): 71-80.
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

