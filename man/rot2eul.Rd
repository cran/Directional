\name{Euler angles  from a rotation matrix on SO(3)}
\alias{rot2eul}
\title{
Compute the Euler angles from a rotation matrix on SO(3).
}

\description{
It calculates three euler angles \eqn{(\theta_{12}, \theta_{13}, \theta_{23})} from a \eqn{(3 \times 3)} rotation matrix X, where X is defined as \eqn{X = R_z(\theta_{12})\times R_y(\theta_{13}) \times R_x(\theta_{23})}. Here \eqn{R_x(\theta_{23})} means a rotation of \eqn{\theta_{23}} radians about the x axis.
}

\usage{
rot2eul(X)
}

\arguments{
\item{X}{
A rotation matrix which is defined as a product of three elementary rotations mentioned above.
Here \eqn{ \theta_{12}, \theta_{23} \in (-\pi, \pi)} and and \eqn{\theta_{13} \in (-\pi/2, \pi/2)}.
}
}

\details{
Given a rotation matrix X, euler angles are computed by equating each element in X with the
corresponding element in the matrix product defined above. This results in nine equations that
can be used to find the euler angles.
}

\value{
For a given rotation matrix, there are two eqivalent sets of euler angles.
}

\references{
Green, P. J. \& Mardia, K. V. (2006). Bayesian alignment using hierarchical models, with applications
in proteins bioinformatics. Biometrika, 93(2):235--254.

http://www.staff.city.ac.uk/~sbbh653/publications/euler.pdf
}

\author{

Anamul Sajib<pmxahsa@nottingham.ac.uk>

R implementation and documentation: Anamul Sajib<pmxahsa@nottingham.ac.uk>
}

\seealso{
\code{\link{eul2rot}
}
}

\examples{

# three euler angles

theta.12 <- sample( seq(-3, 3, 0.3), 1 )
theta.23 <- sample( seq(-3, 3, 0.3), 1 )
theta.13 <- sample( seq(-1.4, 1.4, 0.3), 1 )

theta.12 ; theta.23 ; theta.13

X <- eul2rot(theta.12, theta.23, theta.13)
X  ##  A rotation matrix

e <- rot2eul(X)$v1

theta.12 <- e[3]
theta.23 <- e[2]
theta.13 <- e[1]

theta.12 ; theta.23 ; theta.13

}
