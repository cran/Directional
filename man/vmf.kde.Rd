\name{Mises-Fisher kernel density estimation for (hyper-)spherical data}
\alias{vmf.kde}
\title{
Kernel density estimation for (hyper-)spherical data using a von Mises-Fisher kernel
}

\description{
A von Mises-Fisher kernel is used for the non parametric density estimation.
}

\usage{
vmf.kde(x, h, thumb = "none")
}

\arguments{
\item{x}{
A matrix with unit vectors, i.e. the data being expressed in Euclidean cordinates.
}
  \item{h}{
The bandwidth to be used.
}
  \item{thumb}{
If this is "none", the given bandwidth is used. If it is "rot" the rule of thumb suggested by Garcia-Portugues (2013) is used.
}
}

\details{
A von Mises-Fisher kernel is used for the non parametric density estimation.
}

\value{
A list including:
\item{h}{
The bandwidth used.
}
\item{f}{
A vector with the kernel density estimate calculated for each of the data points.
}
}

\references{
Garcia Portugues, E. (2013). Exact risk improvement of bandwidth selectors for kernel density estimation with directional data. Electronic Journal of Statistics, 7, 1655-1685.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmfkde.tune}, \link{vm.kde}, \link{vmf}, \link{vmkde.tune}
}
}

\examples{
x <- rvmf(100, rnorm(5), 15)
h <- vmfkde.tune(x)[1]
f1 <- vm.kde(x, h = h, thumb = "none")
f2 <- vm.kde(x, h = h, thumb = "rot")
f1$h  ;  f2$h
}

\keyword{ von Mises-Fisher }
\keyword{ Kernel density estimate }
