\name{Contour plot (on the sphere) of the SESPC distribution}
\alias{spher.sespc.contour}
\title{
Contour plot (on the sphere) of the SESPC distribution
}

\description{
The contour plot (on the sphere) of the SESPC distribution is produced.
}

\usage{
spher.sespc.contour(mu, theta, bgcol = "snow", dat = NULL, col = NULL,
lat = 50, long = 50)
}

\arguments{
\item{mu}{
The mean vector the SESPC distribution, a vector in \eqn{R^3}.
}
\item{theta}{
The two \eqn{\theta} parameters of the SESPC distribution.
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
The goal of this function is for the user to see how the SESPC distribution looks like.
}

\value{
A plot containing the contours of the distribution.
}

\references{
Tsagris M. and Alzeley O. (2023). Circular and spherical projected Cauchy distributions: 
A Novel Framework for Circular and Directional Data Modeling.
https://arxiv.org/pdf/2302.02468.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spher.esag.contour}, \link{spher.spcauchy.contour}
}
}

\examples{
\donttest{
mu <- colMeans( as.matrix( iris[, 1:3] ) )
theta <- c(1 ,0.5)
## the lat and long are decreased to 30. Increase them back to 50 to
## see the difference
spher.sespc.contour(mu, theta, lat = 30, long = 30)
}
}

\keyword{ Contour plot}
\keyword{ SESPC distribution }
