\name{Rotation axis and angle of rotation given a rotation matrix}
\alias{Arotation}
\title{
Rotation axis and angle of rotation given a rotation matrix
}

\description{
Given a 3 x 3 rotation matrix, the angle and the axis of rotation are calcualted.
}

\usage{
Arotation(A)
}

\arguments{
  \item{A}{
  A 3 x 3 rotation matrix.
}
}

\details{
If the user does not supply a rotation matrix a message will appear.
}
\value{
A list including:
\item{angle}{
The angle of rotation expressed in degrees.
}
\item{axis}{
The axis of rotation. A vector of two components, latitude and longitude, expressed in degrees.
}
}

\references{
Course webpage of Howard E. Haber.
http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf

Ted Chang (1986). Spherical Regression. Annals of Statistics, 14(3): 907-924.
}

\author{
Michail Tsagris
R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rot.matrix}, \link{rotation}, \link{rsop}
}
}

\examples{
ksi <- c(25.31, 24.29)
theta <- 2.38
A <- rot.matrix(ksi, theta, rads = FALSE)
A
Arotation(A)
}

\keyword{ Angle of rotation }
\keyword{ Axis of rotation }
