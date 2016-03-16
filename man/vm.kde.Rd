\name{vm.kde}
\alias{vm.kde}
\title{
Kernel density estimation of circular data with a von Mises kernel.
}
\description{
Kernel density estimation of circular data with a von Mises kernel.
}
\usage{
vm.kde(u, h, thumb = "none", rads = TRUE)
}
\arguments{
  \item{u}{
  A numeric vector containing the data.
}
  \item{h}{
  The bandwidth.
}
  \item{thumb}{
  It can be either "none", so the bandwidth the user has set will be used, "tay" for the method of Taylor (2008) or "rot" for the method of Garcia-Portugues (2013).
}
  \item{rads}{
  If the data are in radians, this should be TRUE and FALSE otherwise.
}
}
\details{
The user has the option to use a bandwidth he/she has found in some way (cross-validation) or estimate it as Taylor (2008) or Garcia-Portugues (2013).
}
\value{
A list including:
\item{h}{
The bandwidth. If the user chose one of "tay" or "rot" the estimated bandwidth will be returned.
}
\item{f}{
The kernel density estimate at the observed points.
}
}
\references{
Taylor, C. C. (2008). Automatic bandwidth selection for circular density estimation. Computational Statistics & Data Analysis, 52(7): 3493-3500.

Garcia Portugues, E. (2013). Exact risk improvement of bandwidth selectors for kernel density estimation with directional data. Electronic Journal of Statistics, 7, 1655-1685.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou<athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmkde.tune}, \link{vmkde.tune}, \link{vm.kde}
}
}
\examples{
x <- rvonmises(100, 2.4, 10, rads = TRUE)
hist(x, freq = FALSE)
f1 <- vm.kde(x, h = 0.1, thumb = "rot", rads = TRUE)$f
f2 <- vm.kde(x, h = 0.1, thumb = "tay", rads = TRUE)$f
h <- vmkde.tune(x)[1]
f3 <- vm.kde(x, h = h, thumb = "none", rads = TRUE)$f
points(x, f1, col = 1)
points(x, f2, col = 2)
points(x, f3, col = 3)
}
\keyword{ Kernel density}
\keyword{ von Mises kernel }
