\name{Rotation matrix to rotate a spherical vector along the direction of another}
\alias{rotation}
\title{
Rotation matrix to rotate a spherical vector along the direction of another
}

\description{
A rotation matrix is calculated to rotate a unit vector along the direction of another.
}

\usage{
rotation(a, b)
}

\arguments{
\item{a}{
The initial unit vector.
}
\item{b}{
The target unit vector.
}
}

\details{
The function calcualtes a rotation matrix given two vectors. This rotation matrix is the connection between the two spherical only, vectors.
}
\value{
A rotation matrix whose dimension is equal to the length of the unit vectors. 
}

\references{
Amaral G.J.A., Dryden I.L. and Wood A.T.A. (2007).
Pivotal Bootstrap Methods for k-Sample Problems in Directional Statistics and Shape Analysis. 
Journal of the American Statistical Association, 102(478): 695--707.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{Arotation}, \link{rot.matrix}, \link{lambert}, \link{lambert.inv}, \link{rsop}
}
}

\examples{
a <- rnorm(3)
a <- a/sqrt(sum(a^2))
b <- rnorm(3)
b <- b/sqrt(sum(b^2))
A <- rotation(a, b)
A
a  ;  b
a \%*\% t(A)

a <- rnorm(7)
a <- a/sqrt(sum(a^2))
b <- rnorm(7)
b <- b/sqrt(sum(b^2))
A <- rotation(a, b)
A
a  ;  b
a \%*\% t(A)
}
\keyword{ Rotation matrix }
