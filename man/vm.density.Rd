\name{Density of some circular distributions}
\alias{vm.density}
\alias{spml.density}
\alias{wrapcauchy.density}
\alias{circpurka.density}
\title{
Density of some circular distributions
}

\description{
Density of some circular distributions.
}

\usage{
vm.density(x, m, k, rads = FALSE, logden = FALSE)
spml.density(x, mu, rads = FALSE, logden = FALSE)
wrapcauchy.density(x, m, rho, rads = FALSE, logden = FALSE)
circpurka.density(x, m, a, rads = FALSE, logden = FALSE)
}

\arguments{
\item{x}{
A vector with circular data.
}
\item{m}{
The mean value, a scalar. This is the median for the circular Purkayastha distribution.
}
\item{mu}{
The mean vector, a vector with two values.
}
\item{k}{
The concentration parameter.
}
\item{rho}{The \eqn{rho} parameter of the wrapped Cauchy distribution.
}
\item{a}{
The \eqn{alpha} parameter of the circular Purkayastha distribution.
}
\item{rads}{
If the data are in rads, then this should be TRUE, otherwise FALSE.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the von Mises, bivariate projected normal, wrapped Cauchy or the circular Purkayastha distributions is computed.
}

\value{
A vector with the (log) density values of x.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

\seealso{
\code{\link{kent.density}, \link{rvonmises}, \link{ESAGdensity}
}
}

\examples{
x <- rvonmises(500, m = 2.5, k = 10, rads = TRUE)
mod <- circ.summary(x, rads = TRUE, plot = FALSE)
den <- vm.density(x, mod$mesos, mod$kappa, rads = TRUE, logden = TRUE )
mod$loglik
sum(den)
}

