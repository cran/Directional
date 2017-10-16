\name{A test for testing the equality of the concentration parameters for ciruclar data}
\alias{tang.conc}
\title{
A test for testing the equality of the concentration parameter among g samples, where g >= 2 for ciruclar data
}

\description{
A test for testing the equality of the concentration parameter among g samples, where g >= 2 for ciruclar data. It is a tangential approach.
}

\usage{
tang.conc(u, ina, rads = FALSE)
}
\arguments{
\item{u}{
A numeric vector containing the values of all samples.
}
\item{ina}{
A numerical variable or factor indicating the groups of each value.
}
\item{rads}{
If the data are in radians this should be TRUE and FALSE otherwise.
}
}

\details{
This test works for circular data.
}

\value{
A vector including:
\item{test}{
The value of the test statistic.
}
\item{p-value}{
The p-value of the test.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
Fisher, N. I. (1995). Statistical analysis of circular data. Cambridge University Press.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{embed.circaov}, \link{hcf.circaov}, \link{lr.circaov}, \link{het.circaov}, \link{conc.test}
}
}

\examples{
x <- rvonmises(100, 2.4, 15)
ina <- rep(1:4,each = 25)
tang.conc(x, ina, rads = TRUE)
}

\keyword{ Circular data}
\keyword{ Equality of concentrations }
\keyword{ Hypothesis testing }
