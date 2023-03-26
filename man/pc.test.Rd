\name{Hypothesis test for SIPC distribution over the SESPC distribution}
\alias{pc.test}
\title{
Hypothesis test for SIPC distribution over the SESPC distribution
}

\description{
The null hypothesis is whether an SIPC distribution fits the data well, where the altenrative is that SESPC distribution is more suitable.
}

\usage{
pc.test(x, B = 1, tol = 1e-06)
}

\arguments{
\item{x}{
A numeric matrix with three columns containing the data as unit vectors in Euclidean coordinates.
}
\item{B}{
The number of bootstrap re-samples. By default is set to 999. If it is equal to 1, no bootstrap is performed and the
p-value is obtained throught the asymptotic distribution.
}
\item{tol}{
The tolerance to accept that the Newton-Raphson algorithm used in the IAG distribution has converged.
}
}

\details{
Essentially it is a test of rotational symmetry, whether the two \eqn{\theta} parameters are equal to zero.
This works for spherical data only.
}

\value{
A vector including:
\item{test}{
The value of the test statistic.
}
\item{p-value or Bootstrap p-value}{
The p-value of the test.
}
}

\references{
Tsagris M. and Alzeley O. (2023). Circular and spherical projected Cauchy distributions. 
https://arxiv.org/pdf/2302.02468.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{iagesag}, \link{fishkent}, \link{sespc.mle}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
iagesag(x)
pc.test(x)
}

\keyword{ Goodness of fit test }
\keyword{ SESPC distribution }

