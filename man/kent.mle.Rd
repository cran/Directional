\name{MLE of the Kent distribution}
\alias{kent.mle}
\title{
MLe of the Kent distribution
}

\description{
It estimates the concentration and the ovalness parameter of some directional data assuming the Kent distribution. The mean direction and major and minor axes are also estimated.
}

\usage{
kent.mle(x)
}

\arguments{
\item{x}{
A matrix containing spherical data in Euclidean coordinates.
}
}

\details{
The Kent distribution is fitted to some data and its parameters are estimated.
}

\value{
A list including:
\item{runtime}{
The run time of the procedure.
}
\item{G}{
A 3 x 3 matrix whose first column is the mean direction. The second and third columns are the major and minor axes respectively.
}
\item{param}{
A vector with the concentration \eqn{\kappa} and ovalness \eqn{\beta} parameters and the angle \eqn{\psi} used to rotate \bold{H}
and hence estimate \bold{G} as in Kent (1982).
}
\item{logcon}{
The logarithm of the normalising constant, using the third type approximation (Kume and Wood, 2005).
}
\item{loglik}{
The value of the log-likelihood.
}
}

\references{
Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society, Series B, 44(1): 71-80.

Kume Alfred and Wood Andrew T.A. (2005). Saddlepoint approximations for the Bingham and Fisher-Bingham normalizing constants.
Biometrika, 92(2):465-476
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{kent.mle}, \link{fb.saddle}, \link{vmf.mle}, \link{wood.mle}, \link{sphereplot}
}
}

\examples{
x <- rvmf(200, rnorm(3), 15)
kent.mle(x)
vmf.mle(x)
A <- diag( c(-5, 0, 5) )
x <- rfb(200, 15, rnorm(3), A)
kent.mle(x)
vmf.mle(x)
}

\keyword{ Kent distribution }
\keyword{ maximum likelihood estimation }
