\name{vmfkde.tune}
\alias{vmfkde.tune}
\title{
Tuning of the bandwidth parameter in the von Mises-Fisher kernel for (hyper-)spherical data.
}
\description{
Tuning of the bandwidth parameter in the von Mises-Fisher kernel for (hyper-)spherical data. Cross validation is used. 
}
\usage{
vmfkde.tune(x, h = seq(0.1, 1, by = 0.01), plot = TRUE) 
}
\arguments{
  \item{x}{
  A matrix with the data in Euclidean cordinates, i.e. unit vectors.
}
  \item{h}{
  A vector with candidate values for the bandwidth parameter. The search is conducted over this grid of values. 
}
  \item{plot}{
  If TRUE the plot of pseudo-likelihood will appear. 
}
}
\details{
Tuning of the bandwidth parameter in the von Mises-Fisher kernel for (hyper-)spherical data via cross validation. 
}
\value{
A list including:
\item{hopt}{
The optimal bandwidth, the one that maximises the pseudo-likelihood. 
}
\item{cv}{
The value of the maximised pseudo-likelihood. 
}
}
\references{
Garcia Portugues, E. (2013). Exact risk improvement of bandwidth selectors for kernel density estimation with directional data. Electronic Journal of Statistics, 7, 1655--1685.

Wand M. P., & Jones M. C. (1994). Kernel smoothing. CrC Press.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmf.kde}, \link{vmfkde.tune_2},\link{vmf.kerncontour}}
}
\examples{
x <- rvmf(100, rnorm(3), 15)
f1 <- vmfkde.tune(x)
f2 <- vmfkde.tune_2(x)
f1 
f2
}
\keyword{ Kernel density estimate }
\keyword{ von Mises-Fisher distribution }
\keyword{ Tuning of the bandwidth }
