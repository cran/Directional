\name{Tuning of the bandwidth parameter in the von Mises-Fisher kernel}
\alias{vmfkde.tune}
\title{
Tuning of the bandwidth parameter in the von Mises-Fisher kernel for (hyper-)spherical data
}

\description{
Tuning of the bandwidth parameter in the von Mises-Fisher kernel for (hyper-)spherical data whit cross validation.
}

\usage{
vmfkde.tune(x, low = 0.1, up = 1)
}
\arguments{
  \item{x}{
  A matrix with the data in Euclidean cordinates, i.e. unit vectors.
}
  \item{low}{
  The lower value of the bandwdith to search.
}
  \item{up}{
  The upper value of the bandwdith to search.
}
}

\details{
Fast tuning of the bandwidth parameter in the von Mises-Fisher kernel for (hyper-)spherical data via cross validation.
}
\value{
A vector including two elements:
\item{Optimal h}{
The best H found.
}
\item{cv}{
The value of the maximised pseudo-likelihood.
}
}

\references{
Garcia P.E. (2013). Exact risk improvement of bandwidth selectors for kernel density estimation with directional data. Electronic Journal of Statistics, 7, 1655--1685.

Wand M.P. and Jones M.C. (1994). Kernel smoothing. Crc Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.kde},\link{vmf.kerncontour}, \link{vm.kde}, \link{vmkde.tune}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
vmfkde.tune(x)
}

\keyword{ Kernel density estimate }
\keyword{ von Mises-Fisher distribution }
\keyword{ Tuning of the bandwidth }

