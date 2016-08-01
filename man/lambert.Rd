\name{Lambert's equal area projection}
\alias{lambert}
\title{
Lambert's equal area projection
}

\description{
It calculates the Lambert's equal area projection.
}

\usage{
lambert(y)
}

\arguments{
  \item{y}{
  A two column matrix with the data. The first column is the altitude and the second is the longitude.
}
}

\details{
The spherical data are first rotated so that their mean direction is the north pole and then are projectedt on the plane tagent to the sphere at the north pole.
}

\value{
A two-column matrix with the projected points.
}

\references{
Kent, John T. (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society. Series B (Methodological) 44(1):71-80.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{euclid}, \link{lambert.inv}
}
}

\examples{
x <- rvmf(100, rnorm(3), 20)
x <- euclid.inv(x)
a <- lambert(x)
plot(a)
}

\keyword{ Lambert's equal area projection }
\keyword{ Spherical data }
