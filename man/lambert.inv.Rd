\name{Inverse of Lambert's equal area projection}
\alias{lambert.inv}
\title{
Inverse of Lambert's equal area projection
}

\description{
It takes some points from the cartesian coordinates and maps them onto the sphere. The inverse os the Lambert's equal area projection.
}

\usage{
lambert.inv(z, mu)
}

\arguments{
  \item{z}{
  A two- column matrix containing the Lambert's equal rea projected dtaa.
}
  \item{mu}{
The mean direction of the data on the sphere.
}
}

\details{
The data are first mapped on the sphere with mean direction equal to the north pole. Then, they are rotated to have the given mean direction. It is the inverse of the Lambert's equal are projection.
}

\value{
A matrix containing spherical data (unit vectors).
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
\code{\link{lambert}
}
}

\examples{
m <- rnorm(3)
m <- m / sqrt( sum(m^2) )
x <- rvmf(20, m, 19)
mu <- vmf(x)$mu
y <- lambert( euclid.inv(x) )
lambert.inv(y, mu)
euclid.inv(x)
}

\keyword{ Lambert's equal area projection }
\keyword{ Inverse transformation }
