\name{Inverse of the Euclidean transformation}
\alias{euclid.inv}
\title{
Inverse of the Euclidean transformation
}
\description{
It transforms the data from the Euclidan coordinates to latitude dn longitude.
}

\usage{
euclid.inv(U)
}
\arguments{
\item{U}{
A matrix of unit vectors, or even one single unit vector in three dimensions.
}
}

\details{
It takes the matrix of unit vectors and back transforms it to latitude and longitude.
}
\value{
A two column matrix:
\item{u}{ The first column is the latitude and the second is the longitude, both expressed in degrees.
}
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{euclid}, \link{Arotation}, \link{lambert}
}
}

\examples{
x <- rvmf(10, rnorm(3), 10)
euclid.inv(x)
euclid( euclid.inv(x) )
x
}

\keyword{ Euclidean coordinates }
\keyword{ Spherical coordinates }
