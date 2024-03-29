\name{Hypothesis test for IAG distribution over the ESAG distribution}
\alias{iagesag}
\title{
Hypothesis test for IAG distribution over the ESAG distribution
}

\description{
The null hypothesis is whether an IAG distribution fits the data well, where the altenrative is that ESAG distribution is more suitable.
}

\usage{
iagesag(x, B = 1, tol = 1e-07)
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
Essentially it is a test of rotational symmetry, whether the two \eqn{\gamma} parameters are equal to zero.
This works for spherical data only.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
The degrees of freedom of the test. If bootstrap was employed this is "NA".
}
\item{p.value}{
The p-value of the test.
}
\item{alternative}{
A character with the alternative hypothesis.
}
\item{method}{
A character with the test used.
}
\item{data.name}{
A character vector with two elements.
}
}

\references{
Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{fishkent}, \link{iagesag}, \link{pc.test}, \link{esag.mle}, \link{kent.mle},
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
iagesag(x)
fishkent(x, B = 1)
}

\keyword{ Goodness of fit test }
\keyword{ ESAG distribution }
\keyword{ IAG distribution }

