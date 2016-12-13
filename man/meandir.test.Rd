\name{Test for a given mean direction}
\alias{meandir.test}
\title{
Test for a given mean direction
}
\description{
A log-likelihood ratio test for testing whether the sample mena direction is equal to some predefined one.
}
\usage{
meandir.test(x, mu, B = 999)
}
\arguments{
\item{x}{
A matrix with the data, unit vectors.
}
\item{mu}{
A unit vector with the hypothesized mean direction.
}
\item{B}{
A number either 1, so no bootstrap calibration is performed or more than 1, so bootstrap calibration is performed.
}
}
\details{
The log-likelihood ratio test is performed.
}
\value{
A list including:
\item{mean.dir}{
The sample mean direction
}
\item{pvalue}{
The p-value of the test.
}
}
\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmf}, \link{kent.mle}, \link{rayleigh}
}
}
\examples{
mu <- rnorm(5)
mu <- mu / sqrt( sum(mu^2) )
x <- rvmf(100, mu, 10)
meandir.test(x, mu, 1)
meandir.test(x, mu, 499)
}
\keyword{ One sample hypothesis test }
\keyword{ Mean direction }
