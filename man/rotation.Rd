\name{Rotation matrix to rotate a spherical vector along the direction of another}
\alias{rotation}
\title{
Rotation matrix to rotate a spherical vector along the direction of another
}

\description{
A rotation matrix is calcualted to rotate a spherical vector along the direction of another.
}

\usage{
rotation(a, b)
}

\arguments{
\item{a}{
The initial spherical unit vector.
}
\item{b}{
The target spherical unit vector.
}
}

\details{
The function calcualtes a rotation matrix given two vectors. This rotation matrix is the connection between the two spherical only, vectors.
}
\value{
A 3 x 3 matrix with a rotation matrix.
}

\references{
G. J. A. Amaral, I. L. Dryden & Andrew T. A. Wood (2007).
Pivotal Bootstrap Methods for k-Sample Problems in Directional Statistics and Shape Analysis. Journal of the American Statistical Association, 102(478): 695-707.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
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
}
\keyword{ Rotation matrix }
