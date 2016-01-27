\name{mediandir}
\alias{mediandir}
\title{
Spherical and hyperspherical median.
}
\description{
It calculates the (hyper-)spherical median of a sample.
}
\usage{
mediandir(x)
}
\arguments{
  \item{x}{
  The data, a numeric matrix with unit vectors.
}
}
\details{
It employes a numerical optimiser (command optim) to find the median direction as described by Fisher (1985) and Fisher, Lewis and Embleton (1987).
}
\value{
The median direction.
}
\references{
Fisher N. I. (1985). Spherical medians. Journal of the Royal Statistical Society. Series B, 47(2): 342-348.

Fisher N. I., Lewis T. and Embleton B. J. (1987). Statistical analysis of spherical data. Cambridge university press.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmf}, \link{kent.mle}, \link{mediandir_2}
}
}
\examples{
m <- rnorm(3)
m <- m / sqrt( sum(m^2) )
x <- rvmf(100, m, 10)
mediandir(x)
mediandir_2(x)
}
\keyword{ Median direction }
