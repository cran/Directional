\name{euclid}
\alias{euclid}
\title{
Euclidean transformation
}

\description{
It transforms the data from the spherical coordinates to Euclidean coordinates.
}

\usage{
euclid(u)
}

\arguments{
  \item{u}{
  A two column matrix or even one single vector, where the first column (or element) is the latitude and the second is the longitude. The order is important.
}
}

\details{
It takes the matrix of unit vectors of latitude and longitude and  transforms it to unit vectors.
}
\value{
A three column matrix:
\item{U}{ The Euclidean coordinates of the latitude and longitude.
}
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

\seealso{
\code{\link{euclid.inv}, \link{Arotation}, \link{lambert}
}
}

\examples{
x <- rvmf(10, rnorm(3), 10)
u <- euclid.inv(x)
euclid(u)
x
}

\keyword{ Euclidean coordinates }
\keyword{ Spherical coordinates }

