\name{Density of the SESPC distribution}
\alias{dsespc}
\title{
Density of the SESPC distribution
}

\description{
Density of the SESPC distribution.
}

\usage{
dsespc(y, mu, theta, logden = FALSE)
}

\arguments{
\item{y}{
A matrix or a vector with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{mu}{
The mean vector the SESPC distribution, a vector in \eqn{R^3}.
}
\item{theta}{
The two \eqn{\theta} parameters of the SESPC distribution.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the SESPC distribution is computed.
}

\value{
A vector with the (log) density values of y.
}

\references{
Tsagris M. and Alzeley O. (2024). Circular and spherical projected Cauchy distributions: 
A Novel Framework for Circular and Directional Data Modeling.
https://arxiv.org/pdf/2302.02468.pdf

Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{desag}, \link{sespc.mle}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- rsespc(1000, m, c(1, 1))
mod <- sespc.mle(y)
dsespc( y, mod$mu, mod$theta)
}

