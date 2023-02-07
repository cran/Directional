\name{Probability density function of the circular Purkayastha distribution}
\alias{pcircpurka}
\title{
Probability density function of the circular Purkayastha distribution
}

\description{
Probability density function of the circular Purkayastha distribution.
}

\usage{
pcircpurka(u, m, a, rads = FALSE)
}

\arguments{
\item{u}{
A numerical value, either in radians or in degrees.
}
\item{m}{
The mean direction in radians or in degrees.
}
\item{a}{
The \eqn{\alpha} parameter of the circular Purkayastha distribution.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
}

\details{
This value calculates the probability of x being less than theta.
}

\value{
The probability that of x being less than theta, where x follows the circular Purkayastha distribution.
}

\references{
Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through Maximum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70-83

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution. Communications in Statistics-Theory and Methods, 19(6): 1973--1986.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{purka.mle}, \link{dcircpurka}, \link{pvm}
}
}

\examples{
pcircpurka(2, 3, 0.3)
}

