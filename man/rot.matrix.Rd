\name{Rotation matrix from a rotation axis and angle of rotation}
\alias{rot.matrix}
\title{
Rotation matrix from a rotation axis and angle of rotation
}

\description{
It calculates a rotation matrix from a rotation axis and angle of rotation.
}

\usage{
rot.matrix(ksi, theta, rads = FALSE)
}

\arguments{
\item{ksi}{
The rotation axis, a vector with two elements, the first is the
latitude and the second is the longitude.
}
\item{theta}{
The angle of rotation.
}
\item{rads}{
If both the ksi and theta are in rads, this should be TRUE. If both the ksi and theta are in degrees, this should be FALSE.
}
}

\details{
The function accepts as arguments the rotation axis and the angle of rotation and it calcualtes the requested rotation matrix.
}

\value{
A 3 x 3 rotation matrix.
}

\references{
Course webpage of Howard E. Haber.
http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf

Ted Chang (1986). Spherical Regression. Annals of Statistics, 14(3): 907--924.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{Arotation}, \link{rotation}, \link{rsop}
}
}

\examples{
ksi <- c(25.31, 24.29)
theta <- 2.38
A <- rot.matrix(ksi, theta, rads = FALSE)
A
Arotation(A)
}

\keyword{ Rotation matrix }
\keyword{ Axis of rotation }
\keyword{ Angle of rotation }
