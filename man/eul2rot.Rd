\name{Rotation matrix on SO(3) from three Euler angles}
\alias{eul2rot}

\title{
Construct a rotation matrix on SO(3) from the Euler angles.
}

\description{
It forms a rotation matrix X on SO(3) by using three Euler angles \eqn{(\theta_{12}, \theta_{13}, \theta_{23})},
where X is defined as \eqn{X=R_z(\theta_{12}) \times R_y(\theta_{13}) \times R_x( \theta_{23} )}.
Here \eqn{R_x (\theta_{23})} means a rotation of \eqn{\theta_{23}} radians about the x axis.
}

\usage{
eul2rot(theta.12, theta.23, theta.13)
}

\arguments{
\item{theta.12}{
An Euler angle, a number which must lie in \eqn{(-\pi, \pi)}.
}
\item{theta.23}{
An Euler angle, a number which must lie in \eqn{(-\pi, \pi)}.
}
\item{theta.13}{
An Euler angle, a number which must lie in \eqn{(-\pi/2, \pi/2)}.
}
}

\details{
Given three euler angles a rotation matrix X on SO(3) is formed using the transformation according to
Green and Mardia (2006) which is defined above.
}

\value{
A roation matrix.
}

\references{
Green, P. J. \& Mardia, K. V. (2006). Bayesian alignment using hierarchical models, with applications in proteins bioinformatics. Biometrika, 93(2):235--254.

}

\author{

Anamul Sajib<pmxahsa@nottingham.ac.uk>

R implementation and documentation: Anamul Sajib<pmxahsa@nottingham.ac.uk>
}

\seealso{
\code{\link{rot2eul}
}
}

\examples{
# three euler angles

theta.12 <- sample( seq(-3, 3, 0.3), 1 )
theta.23 <- sample( seq(-3, 3, 0.3), 1 )
theta.13 <- sample( seq(-1.4, 1.4, 0.3), 1 )

theta.12 ; theta.23 ; theta.13

X <- eul2rot(theta.12, theta.23, theta.13)
X  # A rotation matrix
det(X)

e <- rot2eul(X)$v1

theta.12 <- e[3]
theta.23 <- e[2]
theta.13 <- e[1]

theta.12 ; theta.23 ; theta.13
}
