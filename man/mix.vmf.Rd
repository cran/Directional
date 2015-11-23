\name{mix.vmf}
\alias{mix.vmf}
\title{
Mixtures of Von Mises-Fisher distributions.
}
\description{
It performs model based clustering for circualr, spherical and hyperspherical data assuming von Mises-Fisher distributions. 
}
\usage{
mix.vmf(x, g)
}
\arguments{
  \item{x}{
  A matrix with the data expressed as unit vectors.
}
  \item{g}{
  The number of groups to fit. It must be greater than or equal to 2. 
}
}
\details{
The initial step of the algorithm is not based on a spherical k-means, but on s imple k-means. The results are comparable to the package movMF. 
}
\value{
A list including:
\item{param}{
A matrix with the mean direction, the concetrations parameter and mixing probability of each group.
}
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{pred}{
The predicted group of each observation. 
}
}
\references{
Kurt Hornik and  Bettina Grun (2014). movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{rmixvmf}, \link{bic.mixvmf}, \link{mixvmf.contour}}
}
\examples{
k <- runif(4, 4, 20)
prob <- c(0.2, 0.4, 0.3, 0.1)
mu <- matrix(rnorm(16), ncol = 4)
mu <- mu / sqrt(rowSums(mu^2))
x <- rmixvmf(200, prob, mu, k)$x
mix.vmf(x, 3)
mix.vmf(x, 4)
mix.vmf(x, 5)
}
\keyword{ Model based clustering}
\keyword{ von Mises-Fisher distribution }
