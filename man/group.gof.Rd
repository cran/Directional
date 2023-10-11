\name{Goodness of fit test for grouped data}
\alias{group.gof}
\title{
Goodness of fit test for grouped data
}

\description{
Goodness of fit test for grouped data.
}

\usage{
group.gof(g, ni, m, k, dist = "vm", rads = FALSE, R = 999, ncores = 1)
}

\arguments{
\item{g}{
A vector with the group points, either in radians or in degrees.
}
\item{ni}{
The frequency of each or group class.
}
\item{m}{
The mean direction in radians or in degrees.
}
\item{k}{
The concentration parameter, \eqn{\kappa}.
}
\item{dist}{
The distribution to be tested, it can be either "vm" or "uniform".
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
\item{R}{
The number of bootstrap simulations to perform, set to 999 by default.
}
\item{ncores}{
The number of cores to use.
}
}

\details{
When you have grouped data, you can test whether the data come from the von Mises-Fisher distribution or from a uniform distribution.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
Since this is a bootstrap based test, there are no degrees of freedom, hence this is "NA".
}
\item{p.value}{
The p-value of the test.
}
\item{alternative}{
A character with the alternative hypothesis.
}
\item{method}{
A character with the test used.
}
\item{data.name}{
A character vector with two elements.
}
}

\references{
Arthur Pewsey, Markus Neuhauser, and Graeme D. Ruxton (2013). Circular Statistics in R.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pvm}, \link{circ.summary}, \link{rvonmises}
}
}

\examples{
x <- rvonmises(100, 2, 10)
g <- seq(min(x) - 0.1, max(x) + 0.1, length = 6)
ni <- as.vector( table( cut(x, g) ) )
group.gof(g, ni, 2, 10, dist = "vm", rads = TRUE, R = 299, ncores = 1)
group.gof(g, ni, 2, 5, dist = "vm", rads = TRUE, R = 299, ncores = 1)
}

