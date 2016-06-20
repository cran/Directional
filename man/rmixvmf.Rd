\name{rmixvmf}
\alias{rmixvmf}
\title{
Random values simulation from a mixture of von Mises-Fisher distributions
}

\description{
The function simulates random values simulation from a given mixture of von Mises-Fisher distributions.
}

\usage{
rmixvmf(n, prob, mu, k)
}

\arguments{
\item{n}{
The sample size.
}
\item{prob}{
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
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mix.vmf}, \link{rvmf}, \link{rvmf}, \link{bic.mixvmf}
}
}

\examples{
k <- runif(4, 4, 20)
prob <- c(0.2, 0.4, 0.3, 0.1)
mu <- matrix(rnorm(16), ncol = 4)
mu <- mu / sqrt( rowSums(mu^2) )
x <- rmixvmf(200, prob, mu, k)$x
bic.mixvmf(x, 5)
}

\keyword{ Random values simulation }
\keyword{ Mixtures of von Mises-Fisher distributions }
