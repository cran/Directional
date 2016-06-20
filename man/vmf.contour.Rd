\name{vmf.contour}
\alias{vmf.contour}
\title{
Contour plots of the von Mises-Fisher distribution on the sphere
}

\description{
Contour plots of the von Mises-Fisher distribution on the sphere.
}

\usage{
vmf.contour(k)
}

\arguments{
\item{k}{
The concentration parameter.
}
}

\details{
The user specifies the concentration parameter only and not the mean direction or data. This is for illustration purposes only. The graph will always contain circles, as the von Mises-Fisher distribution is the analogue of a bivariate normal in two dimensions with a zero covariance.
}

\value{
A contour plot of the von Mises-Fisher distribution.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rvmf}, \link{vmf}, \link{vmf.kerncontour}, \link{kent.contour}
}
}

\examples{
par(mfrow = c(1,3) )
vmf.contour(1)
vmf.contour(5)
vmf.contour(10)
}

\keyword{ von Mises-Fisher distribution }
\keyword{ Contour plot }
