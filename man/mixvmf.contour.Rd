\name{Contour plot of a mixture of von Mises-Fisher distributions model}
\alias{mixvmf.contour}
\title{
Contour plot of a mixture of von Mises-Fisher distributions model for spherical data only.
}
\description{
Contour lines are produced of mixture model for spherical data only.
}
\usage{
mixvmf.contour(u, mod)
}
\arguments{
\item{u}{
A two column matrix. The first column is the longitude and the second is the latitude.
}
\item{mod}{
This is mix.vmf object, actually it is a list. Run a mixture model and save it as mod for example, mod = mix.vmf(x, 3).
}
}
\details{
The contour plot is displayed with latitude and longitude in the axes. No Lambert projection is used here. This works for  spherical data only which are given as longitude and latitude.
}
\value{
A plot including:
The points and the contour lines.
}
\references{
Kurt Hornik and  Bettina Grun (2014). movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmf.kerncontour}, \link{vmf.contour}, \link{mix.vmf}
}
}
\examples{
k <- runif(2, 4, 20)
prob <- c(0.4, 0.6)
mu <- matrix( rnorm(6), ncol = 3 )
mu <- mu / sqrt( rowSums(mu^2) )
x <- rmixvmf(200, prob, mu, k)$x
mod <- mix.vmf(x, 2)
y <- euclid.inv(x)
mixvmf.contour(y, mod)
}
\keyword{ Contour plot }
\keyword{ Mixtures of von Mises-Fisher distributions }
