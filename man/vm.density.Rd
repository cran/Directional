\name{Density of the von Mises distribution}
\alias{vm.density}
\title{
Density of the von Mises distribution
}

\description{
Density of the von Mises distribution.
}

\usage{
vm.density(x, k, m, rads = FALSE, logden = FALSE)
}

\arguments{
\item{x}{
A vector with circular data.
}
\item{k}{
The concentration parameter.
}
\item{m}{
The mean value.
}
\item{rads}{
If the data are in rads, then this should be TRUE, otherwise FALSE.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the von Mises distribution is computed.
}

\value{
A vector with the density values of x.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}
%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{kent.density}, \link{rvonmises}, \link{ESAGdensity}
}
}

\examples{
x <- rvonmises(500, m = 2.5, k = 10, rads = TRUE)
mod <- circ.summary(x, rads = TRUE, plot = FALSE)
den <- vm.density(x, mod$kappa, mod$mesos, rads = TRUE, logden = TRUE )
mod$loglik
sum(den)
}

