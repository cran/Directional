\name{Contour plot of the ESAG distribution without any data}
\alias{esag.contour}
\title{
Contour plot of the ESAG distribution without any data
}

\description{
The contour plot of the ESAG distribution on the sphere is produced. 
}

\usage{
esag.contour(mu, gam, lat, long)
}

\arguments{
\item{mu}{
The mean vector the ESAG distribution, a vector in \eqn{R^3}.
}
\item{gam}{
The two gamma parameters of the ESAG distribution.
}
\item{lat}{
A positive number determing the range of degrees to move left and right from the latitude center. See the example to better understand this argument. 
}
\item{long}{
A positive number determing the range of degrees to move up and down from the longitude center. See the example to better understand this argument.
}
}

\details{
The goal of this function is for the user to see hwo the ESAG distribution looks like.
}

\value{
A plot containing the contours of the ESAG distribution.
}

\references{
Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.
}

\author{
Michail Tsagris and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kent.contour}, \link{vmf.contour}
}
}

\examples{
\donttest{
mu <- colMeans( as.matrix( iris[,1:3] ) )
gam <- c(1,0.5) 
esag.contour(mu, gam, 50, 50)
esag.contour(mu, gam, 30, 40)
}
}

\keyword{ Contour plot}
\keyword{ Kent distribution }
