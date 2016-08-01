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
visual.check(n, F)
}

\arguments{
\item{n}{
The sample size.
}
\item{F}{
An arbitrary 3 x3 matrix represents the parameter matrix of this distribution.
}
}

\details{
For a given parameter matrix F, maximum value of the logarithm of matrix Fisher distribution is calculated via
the form of singular value decomposition of \eqn{F = U \Lambda V^T} which is \eqn{tr(\Lambda)}. Multiply the last
column of \eqn{U} by \eqn{-1} and replace small eigenvalue, say, \eqn{\lambda_3} by \eqn{-\lambda_3} if \eqn{\left| UV^T\right| = -1}.
}

\value{
A plot which shows log probability trace of matrix Fisher distribution.
}

\references{
Habeck M. (2009). Generation of three-dimensional random rotations in fitting and matching problems.
Computational Statistics, 24(4):719--731.
}
\author{

Anamul Sajib<pmxahsa@nottingham.ac.uk>

R implementation and documentation: Anamul Sajib<pmxahsa@nottingham.ac.uk>
}

\examples{
F <- 10^(-1) * matrix( c(85, 11, 41, 78, 39, 60, 43, 64, 48), ncol = 3)

# Here det(UV^{T}) is negative, so small eigenvalue is multiplied by -1.
# The maximum value in this case is (16.17 + 4.80 - 0.57) = 20.4

values <- visual.check(1000, F)

l1 <- values$l1
l2 <- values$l2

l1values <- l1[ seq(1, length(l1), by = 100) ]   # Habeck's method
l2values <- l2[ seq(1, length(l2), by = 100) ]   # Kent's method

plot( 1:length(l1values), l1values, type = "l", ylim = c( min(l1values, l2values ) -0.5,
max(l1values, l2values) + 0.5), col = "red" )
lines( 1:length(l2values), l2values, type = "l", col = "blue" )
abline(h = 20.4, lty = "solid")
}

