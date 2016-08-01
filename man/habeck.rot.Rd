\name{Habeck's rotation matrix generation}
\alias{habeck.rot}
\title{
Generation of three-dimensional random rotations using Habeck's algorithm.
}

\description{
It generates random rotations in three-dimensional space that follow a probability distribution,
matrix Fisher distribution, arising in fitting and matching problem.
}

\usage{
habeck.rot(F)
}

\arguments{
  \item{F}{
  An arbitrary 3 x 3 matrix represents the parameter matrix of this distribution.
}
}

\details{
Firstly rotation matrices \bold{X} are chosen which are the closest to F, and then parameterized using euler angles.
Then a Gibbs sampling algorithm is implemented to generate rotation matrices from the resulting disribution of
the euler angles.
}

\value{
A simulated rotation matrix.
}

\references{
Habeck M (2009). Generation of three-dimensional random rotations in fitting and matching problems. Computational Statistics, 24, 719--731.
}

\author{
Anamul Sajib<pmxahsa@nottingham.ac.uk>

R implementation and documentation: Anamul Sajib<pmxahsa@nottingham.ac.uk>
}

\examples{
F <- 10^(-1) *  matrix( c(85, 11, 41, 78, 39, 60, 43, 64, 48), ncol = 3 )  ## Arbitrary F matrix
X <- habeck.rot(F)
det(X)
}
