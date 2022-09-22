\name{Mixtures of Von Mises-Fisher distributions}
\alias{mixvmf.mle}
\title{
Mixtures of Von Mises-Fisher distributions
}

\description{
It performs model based clustering for circualr, spherical and hyperspherical data assuming von Mises-Fisher distributions.
}

\usage{
mixvmf.mle(x, g, n.start = 10)
}

\arguments{
\item{x}{
A matrix with the data expressed as unit vectors.
}
\item{g}{
The number of groups to fit. It must be greater than or equal to 2.
}
\item{n.start}{
The number of random starts to try. See also R's built-in function \code{\link{kmeans}} for more information about this.
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
\item{iter}{
The number of iteration required by the EM algorithm.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
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
\code{\link{rmixvmf}, \link{bic.mixvmf}, \link{mixvmf.contour}
}
}

\examples{
k <- runif(4, 4, 6)
prob <- c(0.2, 0.4, 0.3, 0.1)
mu <- matrix(rnorm(16), ncol = 4)
mu <- mu / sqrt( rowSums(mu^2) )
x <- rmixvmf(200, prob, mu, k)$x
mixvmf.mle(x, 3)
mixvmf.mle(x, 4)
mixvmf.mle(x, 5)
}

\keyword{ Model based clustering }
\keyword{ von Mises-Fisher distribution }
