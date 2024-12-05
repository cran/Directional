\name{Density of a mixture of rotationally symmetric distributions}
\alias{dmixvmf}
\alias{dmixspcauchy}
\alias{dmixpkbd}
\title{
Density of a mixture of rotationally symmetric distributions
}

\description{
Density of a mixture of rotationally symmetric distributions.
}

\usage{
dmixvmf(y, probs, mu, k, logden = FALSE)
dmixspcauchy(y, probs, mu, rho, logden = FALSE)
dmixpkbd(y, probs, mu, rho, logden = FALSE)
}

\arguments{
\item{y}{
A matrix with unit vectors.
}
\item{probs}{
This is avector with the mixing probability of each group.
}
\item{mu}{
A matrix with the mean direction of each group.
}
\item{k}{
A vector with the concentration parameter of each group.
}
\item{rho}{
A vector with the concentration parameter of each group.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The function computes the density for a given mixture of von Mises-Fisher, spherical Cauchy or Poisson kernel-based distributions.
}

\value{
A vector with the (log) density values of y.
}

\references{
Kurt Hornik and  Bettina Grun (2014). movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf

Tsagris M., Papastamoulis P. and Kato S. (2024). Directional data analysis using the spherical Cauchy and the Poisson kernel-based distribution. https://arxiv.org/pdf/2409.03292.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mixvmf.mle}, \link{rvmf}, \link{bic.mixvmf}
}
}

\examples{
k <- runif(3, 4, 6)
probs <- c(0.2, 0.5, 0.3)
mu <- matrix(rnorm(9), ncol = 3)
mu <- mu / sqrt( rowSums(mu^2) )
x <- rmixvmf(200, probs, mu, k)$x
b <- dmixvmf(x, probs, mu, k)
}
