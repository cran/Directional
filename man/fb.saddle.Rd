\name{Saddlepoint approximations of the Fisher-Bingham distributions}
\alias{fb.saddle}
\title{
Saddlepoint approximations of the Fisher-Bingham distributions
}

\description{
It calculates the logarithm of the normalising constant of the Fisher-Bingham distribution.
}

\usage{
fb.saddle(gam, lam)
}

\arguments{
\item{gam}{
A numeric vector containing the parameters of the Fisher part.
}
\item{lam}{
All the eigenvalues of the Bingham part. Not just the non zero ones.
}
}

\details{
It calculate the three approximations given by Kume and Wood (2005) and it uses the Fisher-Bingham parametrization of that paper.
}

\value{
A list including:
\item{first oder}{
The first order approximation
}
\item{second oder}{
The second order approximation
}
\item{third oder}{
The third order approximation
}
}

\references{
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
\code{\link{kent.logcon}, \link{rfb}, \link{kent.mle}, \link{rbingham}
}
}

\examples{
p <- 3  ;  k <- 1
0.5 * p * log(2 * pi) - (p/2 - 1) * log(k) + log( besselI(k, p/2 - 1, expon.scaled = TRUE) ) + k
## normalising constant of the
## von Mises-Fisher distribution
fb.saddle( c(0, k, 0), c(0, 0, 0) ) ## saddlepoint approximation

## Normalising constant of the Kent distribution
fb.saddle( c(0, 10, 0), c(0, -2, 2) )
kent.logcon(10, 2)
}

\keyword{ Fisher-Bingham distribution }
\keyword{ Normalising constant }
\keyword{ Saddlepoint approximation }

