\name{Contour plot (on the plane) of the ESAG and Kent distributions without any data}
\alias{esag.contour}
\alias{kent.contour}
\title{
Contour plot (on the plane) of the ESAG and Kent and ESAG distributions without any data
}

\description{
The contour plot (on the plane) of the spherical ESAG and Kent distributions is produced.
}

\usage{
esag.contour(mu, gam, lat, long)
kent.contour(k, b)
}

\arguments{
\item{k}{
The concentration parameter.
}
\item{b}{
The ovalness parameter. It has to be less than k/2 in order for the distribution to be unimodal. Otherwise it is bimodal.
}
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
The goal of this function is for the user to see how the Kent or the SAG distribution looks like.
}

\value{
A plot containing the contours of the distribution.
}

\references{
Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society, Series B, 44(1): 71--80.

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
\code{\link{vmf.contour}, \link{vmf.kerncontour}, \link{spher.esag.contour}
}
}

\examples{
\donttest{
kent.contour(10, 4)

mu <- colMeans( as.matrix( iris[,1:3] ) )
gam <- c(1,0.5)
esag.contour(mu, gam, 50, 50)
esag.contour(mu, gam, 30, 40)
}
}

\keyword{ Contour plot}
\keyword{ Kent distribution }
