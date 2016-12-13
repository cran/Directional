\name{MLE of the angular central Gaussian distribution}
\alias{acg}
\title{
MLE of the angular central Gaussian distribution
}

\description{
MLE of the angular central Gaussian distribution
}

\usage{
acg(x)
}

\arguments{
  \item{x}{
  A matrix with directional data, i.e. unit vectors.
}
}

\details{
There is a constraint on the estimated covariance matrix; its trace is equal to the number of variables. An iterative algorithm takes place and convergence is guaranteed.
}

\value{
A list including:
\item{iter}{
The number if iterations required by the algorithm to converge to the solution.
}

\item{cova}{
The estimated covatriance matrix.
}
}

\references{
Tyler D. E. (1987). Statistical analysis for the angular central Gaussian distribution on the sphere.
Biometrika 74(3): 579-589.
}

\author{
Michail Tsagris
R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rbingham}, \link{rfb}, \link{f.rbing}, \link{fb.saddle}
}
}

\examples{
m <- c(0, 0, 0, 0)
s <- cov(iris[, 1:4])
x <- MASS::mvrnorm(500, m, s)
x <- x / sqrt( Rfast::rowsums(x^2) )
mod <- acg(x)
mod
cov2cor(mod$cova)  ## estimated covariance matrix turned into a correlation matrix
cov2cor(s)  ## true covariance matrix turned into a correlation matrix
}

\keyword{ projected normal }
\keyword{ spherical data }
