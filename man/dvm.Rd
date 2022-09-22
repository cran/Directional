\name{Density of some circular distributions}
\alias{dvm}
\alias{dspml}
\alias{dwrapcauchy}
\alias{dcircpurka}
\alias{dggvm}
\alias{dcircbeta}
\alias{dcardio}
\title{
Density of some circular distributions
}

\description{
Density of some circular distributions.
}

\usage{
dvm(x, m, k, rads = FALSE, logden = FALSE)
dspml(x, mu, rads = FALSE, logden = FALSE)
dwrapcauchy(x, m, rho, rads = FALSE, logden = FALSE)
dcircpurka(x, m, a, rads = FALSE, logden = FALSE)
dggvm(x, param, rads = FALSE, logden = FALSE)
dcircbeta(x, m, a, b, rads = FALSE, logden = FALSE)
dcardio(x, m, rho, rads = FALSE, logden = FALSE)
}

\arguments{
\item{x}{
A vector with circular data.
}
\item{m}{
The mean value of the von Mises distribution and of the cardioid, a scalar. This is the median for the circular Purkayastha distribution.
}
\item{mu}{
The mean vector, a vector with two values for the "spml".
}
\item{k}{
The concentration parameter.
}
\item{rho}{The \eqn{rho} parameter of the wrapped Cauchy distribution.
}
\item{a}{
The \eqn{\alpha} parameter of the circular Purkayastha distribution or the \eqn{\alpha} parameter of the
circular beta distribution.
}
\item{b}{
The \eqn{\beta} parameter of the circular beta distribution.
}
\item{param}{
The vector of parameters of the GGVM distribution as returned by the function \code{\link{ggvm.mle}}.
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
\code{\link{dkent}, \link{rvonmises}, \link{desag}
}
}

\examples{
x <- rvonmises(500, m = 2.5, k = 10, rads = TRUE)
mod <- circ.summary(x, rads = TRUE, plot = FALSE)
den <- dvm(x, mod$mesos, mod$kappa, rads = TRUE, logden = TRUE )
mod$loglik
sum(den)
}

