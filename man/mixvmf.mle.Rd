\name{Mixtures of rotationally symmetric distributions}
\alias{mixvmf.mle}
\alias{mixspcauchy.mle}
\alias{mixpkbd.mle}
\title{
Mixtures of rotationally symmetric distributions
}

\description{
It performs model based clustering for circualr, spherical and hyper-spherical data assuming rotationally symetric distributions.
}

\usage{
mixvmf.mle(x, g, n.start = 5, tol = 1e-6, maxiters = 100)
mixspcauchy.mle(x, g, n.start = 5, tol = 1e-6, maxiters = 100)
mixpkbd.mle(x, g, n.start = 5, tol = 1e-6, maxiters = 100)
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
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxiters}{
The maximum number of iterations to perform.
}
}

\details{
The initial step of the algorithm is not based on a spherical k-means, but on simple k-means. The results are comparable to the package movMF for the mixtures of von Mises-Fisher distributions. The other cases are mixtures of spherical Cauchy distributions or mixtures of Poisson kernel-based distributions.
}

\value{
A list including:
\item{param}{
A matrix with the mean direction, the concentration parameters and mixing probability of each group.
}
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{pred}{
The predicted group of each observation.
}
\item{w}{
The estimated probabilities of each observation to belong to each cluster.
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

Tsagris M., Papastamoulis P. and Kato S. (2024). Directional data analysis using the spherical Cauchy and the Poisson kernel-based distribution. https://arxiv.org/pdf/2409.03292.
}

\author{
Michail Tsagris and Panagiotis Papastamoulis.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Panagiotis Papastamoulis \email{papastamoulis@aueb.gr}.
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
