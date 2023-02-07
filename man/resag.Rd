\name{Simulation of random values from the ESAG distribution}
\alias{resag}
\title{
Simulation of random values from the ESAG distribution
}

\description{
Simulation of random values from the ESAG distribution.
}

\usage{
resag(n, mu, gam)
}

\arguments{
\item{n}{
A number; how many vectors you want to generate.
}
\item{mu}{
The mean vector the ESAG distribution, a vector in \eqn{R^3}.
}
\item{gam}{
The two \eqn{\gamma} parameters of the ESAG distribution.
}
}

\details{
A random sample from the ESAG distribution is generated. In case the \eqn{\gamma_s} are zero, the sample is drawn
from the Independent Angular Gaussian (or projected normal).
}

\value{
An \eqn{n \times 3} matrix with the simulated unit vectors.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{esag.mle}, \link{desag}, \link{spml.mle}, \link{acg.mle}, \link{circ.summary}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- resag(1000, m, c(1, 0.5) )
esag.mle(y)
}

\keyword{ ESAG distribution }
\keyword{ simulation }
\keyword{ spherical data }
