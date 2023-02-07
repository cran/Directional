\name{Simulation of random values from a mixture of von Mises-Fisher distributions}
\alias{rmixvmf}
\title{
Simulation of random values from a mixture of von Mises-Fisher distributions
}

\description{
The function simulates random values simulation from a given mixture of von Mises-Fisher distributions.
}

\usage{
rmixvmf(n, probs, mu, k)
}

\arguments{
\item{n}{
The sample size.
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
}

\details{
The function simulates random values simulation from a given mixture of von Mises-Fisher distributions using the rvmf function.
}

\value{
A list including:
\item{id}{
An indicator of the group of each simulated vector.
}
\item{x}{
A matrix with the simulated data.
}
}

\references{
Kurt Hornik and  Bettina Grun (2014). movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf}

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
k <- runif(3, 4, 20)
probs <- c(0.2, 0.5, 0.3)
mu <- matrix(rnorm(9), ncol = 3)
mu <- mu / sqrt( rowSums(mu^2) )
x <- rmixvmf(200, probs, mu, k)$x
bic.mixvmf(x, 5)
}

\keyword{ Random values simulation }
\keyword{ Mixtures of von Mises-Fisher distributions }
