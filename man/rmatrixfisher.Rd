\name{Simulation from a Matrix Fisher distribution on SO(3)}
\alias{rmatrixfisher}
\title{
Simulation from a Matrix Fisher distribution on SO(3)
}

\description{
It simulates random samples (rotation matrices) from a Matrix Fisher distribution with any given parameter matrix, F (3x3).
}

\usage{
rmatrixfisher(n, F)
}

\arguments{
\item{n}{
the sample size.
}
\item{F}{
An arbitrary 3x3 matrix.
}
}

\details{
Firstly corresponding Bingham parameter A is determined for a given Matrix Fisher parameter F using John Kent (2013) algorithm and then Bingham samples for parameter A are generated using rbingham code. Finally convert Bingham samples to Matrix Fisher samples according to the Kent (2013) transformation.
}

\value{
An array with simulated rotation matrices.
}

\references{
Kent J.T., Ganeiber A.M. and Mardia K.V. (2013). A new method to simulate the Bingham and related
distributions in directional data analysis with applications. http://arxiv.org/pdf/1310.8110v1.pdf
}

\author{
Anamul Sajib & Chris Fallaize.

R implementation and documentation: Anamul Sajib <sajibstat@du.ac.bd> & Chris Fallaize.
}

\seealso{
\code{\link{matrixfisher.mle}
}
}

\examples{
F <- matrix( c(85, 11, 41, 78, 39, 60, 43, 64, 48), ncol = 3) / 10   ### An arbitrary F matrix
a <- rmatrixfisher(10, F)
}
