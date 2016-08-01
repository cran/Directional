\name{Contour plot of spherical data using a von Mises-Fisher kernel density estimate}
\alias{vmf.kerncontour}
\title{
Contour plot of spherical data using a von Mises-Fisher kernel density estimate
}

\description{
Contour plot of spherical data using a von Mises-Fisher kernel density estimate.
}

\usage{
vmf.kerncontour(u, thumb = "none")
}

\arguments{
\item{u}{
A two column matrix. The first coolumn is the latitude and the second is the longitude.
}
\item{thumb}{
This is either 'none' (defualt), or 'rot' for the rule of thumb suggested by Garcia-Portugues (2013). If it is "none" it is estimated via cross validation, with the fast function "vmfkde.tune_2".
}
}

\details{
It calcculates the contour plot using a von Mises-Fisher kernel for spherical data only.
}

\value{
A plot with the data in latitude and longitude along with their kernel density estimate contour lines.
}

\references{
Garcia Portugues, E. (2013). Exact risk improvement of bandwidth selectors for kernel density estimation with directional data. Electronic Journal of Statistics, 7, 1655--1685.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.kde}, \link{vmfkde.tune}, \link{vmf.contour}, \link{kent.datacontour}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
x <- euclid.inv(x)
par( mfrow = c(1, 2) )
vmf.kerncontour(x, "rot")
vmf.kerncontour(x, "none")
}

\keyword{ Contour plot }
\keyword{ von Mises-Fisher kernel }
