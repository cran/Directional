\name{Converting an unsigned unit quaternion to rotation matrix on SO(3)}
\alias{quat2rot}
\title{
Converting an unsigned unit quaternion to rotation matrix on SO(3)
}

\description{
It forms a (3 x 3) rotation matrix on SO(3) from an unsigned unite quaternion in \eqn{S^3} (the four-dimensional sphere).
}

\usage{
quat2rot(x)
}

\arguments{
\item{x}{
An unsigned unit quaternion in \eqn{S^3}.
}
}

\details{
Given an unsigned unit quaternion in \eqn{S^3} it forms a rotation matrix on SO(3), according to the transformation proposed by Prentice (1986).
}

\value{
A rotation matrix.
}

\references{
Prentice,M. J. (1986). Orientation statistics without parametric assumptions.Journal of the Royal Statistical Society. Series B: Methodological 48(2).
}

\author{
Anamul Sajib

R implementation and documentation: Anamul Sajib <sajibstat@du.ac.bd>
}

\seealso{
\code{\link{rot2quat}, \link{rotation}, \link{Arotation} \link{rot.matrix}
}
}

\examples{
x <- rnorm(4)
x <- x/sqrt( sum(x^2) )
x                   ## an unit quaternion in R4 ##
quat2rot(x)
}
