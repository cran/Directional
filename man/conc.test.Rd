\name{Test of equality of the concentration parameters for circular data}
\alias{conc.test}
\title{
A test for testing the equality of the concentration parameter among g samples, where g >= 2 for ciruclar data
}

\description{
A test for testing the equality of the concentration parameter among g samples, where g >= 2 for ciruclar data.
}

\usage{
conc.test(u, ina, rads = FALSE)
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
A list including:
\item{mess}{
A message informing the use of the test statistic used.
}
\item{res}{
A numeric vector containing the value of the test statistic and its associated p-value.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}
\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{embed.circaov}, \link{hcf.circaov}, \link{lr.circaov}, \link{het.circaov}
}
}

\examples{
x <- rvonmises(100, 2.4, 15)
ina <- rep(1:4,each = 25)
conc.test(x, ina, rads = TRUE)
}

\keyword{ Circular data }
\keyword{ Equality of concentrations }
\keyword{ Hypothesis testing }
