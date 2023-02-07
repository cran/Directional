\name{Simulation from a Bingham distribution using any symmetric matrix A}
\alias{rbingham}
\title{
Simulation from a Bingham distribution using any symmetric matrix A
}

\description{
It simulates random values from a Bingham distribution with any given symmetric matrix.
}
\usage{
rbingham(n, A)
}
\arguments{
  \item{n}{
  The sample size.
}
  \item{A}{
  A symmetric matrix.
}
}

\details{
The eigenvalues are fist calcualted and then Chris Fallaize and Theo Kypraio's code (f.rbing) is used. The resulting simulated data anre then right multiplied by the eigenvectors of the matrix A.
}

\value{
A matrix with the siumlated data.
}

\references{
Kent J. T., Ganeiber A. M. & Mardia K. V. (2018). A new unified approach for the simulation of a wide class of directional distributions. 
Journal of Computational and Graphical Statistics, 27(2): 291--301.

Kent J.T., Ganeiber A.M. and Mardia K.V. (2013).
A new method to simulate the Bingham and related distributions
in directional data analysis with applications.
http://arxiv.org/pdf/1310.8110v1.pdf

Fallaize C. J. & Kypraios T. (2016). Exact bayesian inference for the Bingham distribution. 
Statistics and Computing, 26(1): 349--360.
http://arxiv.org/pdf/1401.2894v1.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{f.rbing}, \link{rfb}, \link{rvmf}, \link{rkent}
}
}

\examples{
A <- cov(iris[, 1:3])
x <- rbingham(100, A)
}

\keyword{ Bingham distribution }
\keyword{ Simulation of random values }
