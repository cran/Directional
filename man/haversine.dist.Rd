\name{Haversine distance matrix}
\alias{haversine.dist}
\title{
Harvesine distance matrix
}

\description{
Haversine distance matrix.
}

\usage{
haversine.dist(x)
}

\arguments{
\item{x}{
A a matrix of two columns. The first column is the latitude and the second the longitude.
}
}

\details{
The function computes the haversine distance between all observations.
}

\value{
A matrix with the haversine distances between all observations.
}

\references{
https://en.wikipedia.org/wiki/Haversine_formula
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cosnn}, \link{dirknn}
}
}

\examples{
x <- rvmf(10, rnorm(3), 10)
x <- euclid.inv(x)
haversine.dist(x)
}
