\name{vmfkde.tune_2}
\alias{vmfkde.tune_2}
\title{
Fast tuning of the bandwidth parameter in the von Mises-Fisher kernel for (hyper-)spherical data with an optimizer.
}
\description{
Fast tuning of the bandwidth parameter in the von Mises-Fisher kernel for (hyper-)spherical data when cross validation is used with an optimiser. 
}
\usage{
vmfkde.tune_2(x, low = 0.1, up = 1)
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
Garcia Portugues, E. (2013). Exact risk improvement of bandwidth selectors for kernel density estimation with directional data. Electronic Journal of Statistics, 7, 1655--1685.

Wand, M. P., & Jones, M. C. (1994). Kernel smoothing. Crc Press.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmf.kde}, \link{vmfkde.tune},\link{vmf.kerncontour}}
}
\examples{
x <- rvmf(100, rnorm(3), 15)
f1 <- vmfkde.tune(x)
f2 <- vmfkde.tune_2(x)
f1  ;  f2
}
\keyword{ Kernel density estimate }
\keyword{ von Mises-Fisher distribution }
\keyword{ Tuning of the bandwidth }

