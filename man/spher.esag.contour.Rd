\name{Contour plot (on the sphere) of the ESAG and Kent distributions}
\alias{spher.esag.contour}
\alias{spher.kent.contour}
\title{
Contour plot (on the sphere) of the ESAG and Kent distributions
}

\description{
The contour plot (on the sphere) of the ESAG and Kent distributions is produced.
}

\usage{
spher.esag.contour(mu, gam, bgcol = "snow", dat = NULL, col = NULL,
lat = 50, long = 50)
spher.kent.contour(G, param, bgcol = "snow", dat = NULL, col = NULL,
lat = 50, long = 50)
}

\arguments{
\item{mu}{
The mean vector the ESAG distribution, a vector in \eqn{R^3}.
}
\item{gam}{
The two gamma parameters of the ESAG distribution.
}
\item{G}{
For the Kent distribution, a 3 x 3 matrix whose first column is the mean direction. The second and third columns are the major and minor axes respectively.
}
\item{param}{
For the Kent distribution a vector with the concentration \eqn{\kappa} and ovalness \eqn{\beta} parameters. The angle \eqn{\psi} has been absorbed inside the matrix G.
}
\item{bgcol}{
The color of the surface of the sphere.
}
\item{dat}{
If you have you want to plot supply them here. This has to be a numerical matrix with three columns, i.e. unit vectors.
}
\item{col}{
If you supplied data then choose the color of the points. If you did not choose a color, the points will appear in red.
}
\item{lat}{
A positive number determing the range of degrees to move left and right from the latitude center. See the example to better understand this argument.
}
\item{long}{
A positive number determing the range of degrees to move up and down from the longitude center. See the example to better understand this argument.
}
}

\details{
The goal of this function is for the user to see how the ESAG or the Kent distribution looks like.
}

\value{
A plot containing the contours of the distribution.
}

\references{
Kent John (1982). The Fisher-Bingham distribution on the sphere.
Journal of the Royal Statistical Society, Series B, 44(1): 71--80.

Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{esag.contour}, \link{spher.purka.contour}, \link{kent.contour}
}
}

\examples{
\donttest{
mu <- colMeans( as.matrix( iris[, 1:3] ) )
gam <- c(1 ,0.5)
## the lat and long are decreased to 30. Increase them back to 50 to
## see the difference
spher.esag.contour(mu, gam, lat = 30, long = 30)
}
}

\keyword{ Contour plot}
\keyword{ ESAG distribution }
\keyword{ Kent distribution }
