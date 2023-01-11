\name{Converting a rotation matrix on SO(3) to an unsigned unit quaternion}
\alias{rot2quat}
\title{
Converting a rotation matrix on SO(3) to an unsigned unit quaternion
}

\description{
It returns an unsigned unite quaternion in \eqn{S^3} (the four-dimensional sphere) from a \eqn{3 \times 3} 
rotation matrix on SO(3).
}

\usage{
rot2quat(X)
}

\arguments{
\item{X}{
A rotation matrix in SO(3).
}
}

\details{
Firstly construct a system of linear equations by equating the corresponding components of the theoretical rotation 
matrix proposed by Prentice (1986), and given a rotation matrix. Finally, the system of linear equations are solved 
by following the tricks mentioned in second reference here in order to achieve numerical accuracy to get quaternion values.
}

\value{
A unsigned unite quaternion.
}

\references{
Prentice,M. J. (1986). Orientation statistics without parametric assumptions.Journal of the 
Royal Statistical Society. Series B: Methodological 48(2).  //http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
}

\author{
Anamul Sajib.

R implementation and documentation: Anamul Sajib <sajibstat@du.ac.bd>.
}

\seealso{
\code{\link{quat2rot}, \link{rotation}, \link{Arotation} \ link{rot.matrix}
}
}

\examples{
x <- rnorm(4)
x <- x/sqrt( sum(x^2) ) ## an unit quaternion in R4 ##
R <- quat2rot(x)
R
x
rot2quat(R) ## sign is not exact as you can see
}
