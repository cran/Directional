\name{Probability density function of the von Mises-Fisher distribution}
\alias{pvm}
\title{
Probability density function of the von Mises-Fisher distribution
}

\description{
Probability density function of the von Mises-Fisher distribution.
}

\usage{
pvm(theta, m, k, rads = FALSE)
}

\arguments{
\item{theta}{
A numerical value, either in radians or in degrees.
}
\item{m}{
The mean direction in radians or in degrees.
}
\item{k}{
The concentration parameter, \eqn{\kappa}.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
}

\details{
This value calculates the probability of x being less than theta and is used by \code{\link{group.gof}}.
}

\value{
The probability that of x being less than theta, where x follows the von Mises-Fisher distribution.
}

\references{
Arthur Pewsey, Markus Neuhauser, and Graeme D. Ruxton (2013).
Circular Statistics in R.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{group.gof}, \link{circ.summary}, \link{rvonmises}
}
}

\examples{
pvm(1, 2, 10, rads = TRUE)
pvm(2, 2, 10, rads = TRUE)
}

\keyword{ Probability density function }
\keyword{ von Mises-Fisher distribution }
