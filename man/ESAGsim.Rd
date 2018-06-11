\name{Random values generation from the ESAG distribution}
\alias{ESAGsim}
\title{
Random values generation from the ESAG distribution
}

\description{
Random values generation from the ESAG distribution.
}

\usage{
ESAGsim(n, param)
}

\arguments{
\item{n}{
A number; how many vectors you want to generate.
}
\item{param}{
The parameters of the ESAG distribution, the first three are the mean vector in \eqn{R^3} and the next
two are the two gammas.
}
}

\details{
A random sample from the ESAG distribution is generated. In case the gammas are zeo the sample is drawn
from Independent Angular Gaussian (IAG) or projected normal.
}

\value{
An nx3 matrix with the simualted unit vectors.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

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
\code{ \link{ESAGmle}, \link{ESAGdensity}, \link{spml.mle}, \link{acg}, \link{circ.summary}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- ESAGsim(1000, c(m, 1, 0.5) )
ESAGmle(y)
}

\keyword{ ESAG distribution }
\keyword{ simulation }
\keyword{ spherical data }
