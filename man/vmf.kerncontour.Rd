\name{Contour plot of spherical data using a von Mises-Fisher kernel density estimate}
\alias{vmf.kerncontour}
\title{
Contour plot of spherical data using a von Mises-Fisher kernel density estimate
}

\description{
Contour plot of spherical data using a von Mises-Fisher kernel density estimate.
}

\usage{
vmf.kerncontour(u, thumb = "none", den.ret = FALSE, full = FALSE, ngrid = 100)
}

\arguments{
\item{u}{
A two column matrix. The first coolumn is the latitude and the second is the longitude.
}
\item{thumb}{
This is either 'none' (defualt), or 'rot' for the rule of thumb suggested by Garcia-Portugues (2013).
If it is "none" it is estimated via cross validation, with the fast function \code{\link{vmfkde.tune}}.
}
\item{den.ret}{
If FALSE (default), plots the contours of the density along with the individual
points. If TRUE, will instead return a list with the Longitudes, Latitudes and
Densities. Look at the 'value' section for details.
}
\item{full}{
If FALSE (default), uses the range of positions from 'u' to calculate and
optionally plot densities. If TRUE, calculates densities covering the entire
sphere.
}
\item{ngrid}{
Sets the resolution of the density calculation.
}
}

\details{
It calculates the contour plot using a von Mises-Fisher kernel for spherical data only.
}

\value{
The contour lines of the data. If "den.ret" was set to TRUE a list including:
\item{lat}{
The latitude values.
}
\item{long}{
The longitude values.
}
\item{h}{
The optimal bandwidth.
}
\item{den}{
The kernel density estimate contour points.
}
}

\references{
Garcia Portugues, E. (2013). Exact risk improvement of bandwidth selectors for kernel density estimation with directional data.
Electronic Journal of Statistics, 7, 1655--1685.
}

\author{
Michail Tsagris, Micah J. Waldstein and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr},
Micah J. Waldstein \email{micah@waldste.in} and Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.kde}, \link{vmfkde.tune}, \link{vmf.contour}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
x <- euclid.inv(x)
\donttest{
vmf.kerncontour(x, "rot")
}
}

\keyword{ Contour plot }
\keyword{ von Mises-Fisher kernel }
