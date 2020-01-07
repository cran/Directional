\name{Check visually whether matrix Fisher samples is correctly generated or not}
\alias{visual.check}
\title{
Check visually whether matrix Fisher samples is correctly generated or not.
}

\description{
It plots the log probability trace of matrix Fisher distribution which should close to the maximum value of the
logarithm of matrix Fisher distribution, if samples are correctly generated.
}

\usage{
visual.check(x, Fa)
}

\arguments{
\item{x}{
The simulated data. An array with at least 2 3x3 matrices.
}
\item{Fa}{
An arbitrary 3x3 matrix represents the parameter matrix of this distribution.
}
}

\details{
For a given parameter matrix Fa, maximum value of the logarithm of matrix Fisher distribution is calculated via
the form of singular value decomposition of \eqn{Fa = U \Lambda V^T} which is \eqn{tr(\Lambda)}. Multiply the last
column of \eqn{U} by \eqn{-1} and replace small eigenvalue, say, \eqn{\lambda_3} by \eqn{-\lambda_3} if \eqn{| UV^T| = -1}.
}

\value{
A plot which shows log probability trace of matrix Fisher distribution. The values are also returned.
}

\references{
Habeck M. (2009). Generation of three-dimensional random rotations in fitting and matching problems.
Computational Statistics, 24(4):719--731.
}
\author{

Anamul Sajib<sajibstat@du.ac.bd>

R implementation and documentation: Anamul Sajib<sajibstat@du.ac.bd>
}

\examples{
Fa <- matrix( c(85, 11, 41, 78, 39, 60, 43, 64, 48), ncol = 3) / 10
x <- rmatrixfisher(1000, Fa)
a <- visual.check(x, Fa)
}

