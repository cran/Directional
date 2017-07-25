\name{Density of the ESAG distribution}
\alias{ESAGdensity}
\title{
Density of the ESAG distribution
}

\description{
Density of the ESAG distribution.
}

\usage{
ESAGdensity(y, param, logden = FALSE)
}

\arguments{
\item{y}{
A matrix or a vector with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{param}{
The parameters of the ESAG distribution, the first three are the mean vector in \eqn{R^3} and the next two are the two gammas.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
MLE of the projected normal distribution, on the sphere, is implemented.
}

\value{
A vector with the density values of y.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Paine P.J., Preston S.P., Tsagris M and Wood A.T.A. (2017). An Elliptically Symmetric Angular Gaussian Distribution. Statistics and Computing (Accepted).
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ESAGmle}, \link{ESAGsim}, \link{spml.mle}, \link{acg}, \link{circ.summary}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- ESAGsim(1000, c(m, 1,0.5) )
mod <- ESAGmle(y)
ESAGdensity( y, c(mod$mu, mod$gam) )
}

\keyword{ ESAG distribution }
\keyword{ spherical data }
