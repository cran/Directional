\name{Simulation of random values from the SESPC distribution}
\alias{rsespc}
\title{
Simulation of random values from the SESPC distribution
}

\description{
Simulation of random values from the SESPC distribution
}

\usage{
rsespc(n, mu, theta)
}

\arguments{
\item{n}{
A number; how many vectors you want to generate.
}
\item{mu}{
The mean vector the SESPC distribution, a vector in \eqn{R^3}.
}
\item{theta}{
The two \eqn{\theta} parameters of the SESPC distribution.
}
}

\details{
A random sample from the SESPC distribution is generated. In case the \eqn{\theta_s} are zero, the sample is drawn
from the SIPC (spherical independent projected Cauchy) distribution.
}

\value{
An \eqn{n \times 3} matrix with the simulated unit vectors.
}

\references{

Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{sespc.mle}, \link{dsespc}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- rsespc(1000, m, c(1, 0.5) )
sespc.mle(y)
}

\keyword{ ESAG distribution }
\keyword{ simulation }
\keyword{ spherical data }
