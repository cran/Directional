\name{MLE of the Wood bimodal distribution on the sphere}
\alias{wood.mle}
\title{
MLE of the Wood bimodal distribution on the sphere
}

\description{
It estimates the parameters of the Wood bimodal distribution.
}

\usage{
wood.mle(y)
}

\arguments{
\item{y}{
A matrix containing two columns. The first one is the latitude and the second is the longitude, both expressed in degrees.
}
}

\details{
The Wood distribution is fitted to some data and its parameters are estimated. It is a bimodal distribution which contains 5 parameters, just like the Kent distribution.
}

\value{
A list including:
\item{info}{
A 5 x 3 matrix containing the 5 parameters, \eqn{\gamma}, \eqn{\delta}, \eqn{\alpha}, \eqn{\beta} and \eqn{\kappa} along with their corresponding 95\% confidence intervals all expressed in degrees.
}
\item{modes}{
The two axis of the modes of the distribution expressed in degrees.
}
\item{unitvectors}{
A 3 x 3 matrix with the 3 unit vectors associated with the \eqn{\gamma} and \eqn{\delta} parameters.
}
\item{loglik}{
The value of the log-likelihood.
}
}

\references{
Wood A.T.A. (1982). A bimodal distribution on the sphere.
Journal of the Royal Statistical Society, Series C, 31(1): 52--58.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kent.mle}, \link{esag.mle}, \link{vmf.mle}, \link{sphereplot}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
x <- euclid.inv(x)
wood.mle(x)
}

\keyword{ Wood distribution }
\keyword{ Bimodal distribution on the sphere }
