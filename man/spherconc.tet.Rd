\name{spherconc.test}
\alias{spherconc.test}
\title{
Test for equality of concentration parameters for spherical data
}
\description{
This tests the equality of concentration parameters for spherical data only. 
}
\usage{
spherconc.test(x, ina)
}
\arguments{
  \item{x}{
  A matrix with the data in Euclidean coordinates, i.e. unit vectors
}
  \item{ina}{
  A variable indicating the groupings of the observations.
}
}
\details{
The test is designed for spherical data only. 
}
\value{
A list including:
\item{mess}{
A message stating the value of the mean resultant and which test statistic was used, U1, U2 or U3. 
}
\item{res}{
A vector containing the test statistic and its p-value.
}
}
\references{
Kanti V. Mardia and Peter E. Jupp. Directional statistics, p.g. 226-227. 
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{het.aov}, \link{lr.aov}, \link{embed.aov}, \link{hcf.aov} }
}
\examples{
x <- rvmf(100, rnorm(3), 15)
ina <- rep(1:4, each = 25)
spherconc.test(x, ina)
}
\keyword{ Concentration parameters }
\keyword{ Hypothesis testing }
\keyword{ Spherical data }

