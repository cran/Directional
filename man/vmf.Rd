\name{MLE of the von Mises-Fisher distribution}
\alias{vmf}
\title{
MLE of von Mises-Fisher distribution
}

\description{
MLE of the von Mises-Fisher distribution.
}

\usage{
vmf(x, fast = FALSE, tol = 1e-07)
}

\arguments{
\item{x}{
A matrix with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{fast}{
A boolean variable to do a faster implementation.
}
\item{tol}{
The tolerance to accept that the E-M algorithm used to estimate the concentration parameter has converged.
}
}

\details{
The mean direction and concentration of a fitted von Mises-Fisher distribution are estimated.
}

\value{
If fast = FALSE a list including all the following. If fast = TRUE less items are returned.
\item{mu}{
The mean direction.
}
\item{kappa}{
The concentration parameter.
}
\item{MRL}{
The mean resultant length.
}
\item{vark}{
The variance of the concentration parameter.
}
\item{loglik}{
The maximum log-likelihood value.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra, S. (2012). A short note on parameter approximation for von Mises-Fisher distributions: and a fast implementation of Is(x). Computational Statistics, 27(1): 177--190.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{iag.mle}, \link{rvmf}, \link{kent.mle}, \link{vmf.kde}, \link{wood.mle}
}
}

\examples{
m <- rnorm(3)
m <- m/sqrt(sum(m^2))
m
x <- rvmf(100, m, 7)
vmf(x)
x <- rvmf(500, m, 7)
vmf(x)
}

\keyword{ Maximum likelihood estimation }
\keyword{ von Mises-Fisher distribution }
