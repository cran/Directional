\name{vmkde.tune}
\alias{vmkde.tune}
\title{
Tuning of the bandwidth parameter in the von Mises kernel for circular data.
}
\description{
Tuning of the bandwidth parameter in the von Mises kernel for circular data. Cross validation is used.
}
\usage{
vmkde.tune(u, h = seq(0.1, 1, by = 0.01), rads =T, plot = TRUE, ncores = 4)
}
\arguments{
  \item{u}{
  The data, a numerical vector.
}
  \item{h}{
    A vector with candidate values for the bandwidth parameter. The search is conducted over this grid of values.
}
  \item{rads}{
  If the data are in radians this should be TRUE and FALSE otherwise.
}
  \item{plot}{
  If TRUE the plot of pseudo-likelihood will appear.
}
 \item{ncores}{
  How many cores should be used? If more than 1, the computations take place in parallel.
}
}
\details{
Tuning of the bandwidth parameter in the von Mises kernel for circula data via cross validation.
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
Taylor C. C. (2008). Automatic bandwidth selection for circular density estimation. Computational Statistics & Data Analysis, 52(7), 3493--3500.

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
\code{\link{vm.kde}, \link{vmkde.tune_2}
}
}
\examples{
u <- rvonmises(100, 2.4, 10, rads = TRUE)
vmkde.tune(u, ncores = 2)
vmkde.tune_2(u)
}
\keyword{ Kernel density estimate }
\keyword{ von Mises distribution }
\keyword{ Tuning the bandwidth }
