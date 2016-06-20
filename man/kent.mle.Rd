\name{kent.mle}
\alias{kent.mle}
\title{
MLe of the Kent distribution.
}

\description{
It estimates the concentration and the ovalness parameter of some directional data assuming the Kent distribution. The mean direction and major and minor axes are also estimated.
}

\usage{
kent.mle(x)
}

\arguments{
  \item{x}{
A matrix containing spherical data; in Euclidean coordinates.
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
\item{para}{
The concentration \eqn{\kappa} and ovalness \eqn{\beta} parameters are stored here.
%,along with the angle \eqn{\psi} used to rotate \bold{H} and hence estimate \eqn{\textbf{\Gamma}} in Kent (1982).
}
\item{logcon}{
The logarithm of the normalising constant, using the third type apprxoimation.
}
\item{loglik}{
The value of the log-likelihood.
}
}

\references{
Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society, Series B, 44(1): 71-80.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{kent.mle}, \link{fb.saddle}, \link{vmf}, \link{wood.mle}
}
}

\examples{
x <- rvmf(200, rnorm(3), 15)
kent.mle(x)
vmf(x)
A <- diag( c(-5, 0, 5) )
x <- rfb(200, 15, rnorm(3), A)
kent.mle(x)
vmf(x)
}

\keyword{ Kent distribution }
\keyword{ Maximum likelihood estimation }
