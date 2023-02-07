\name{Hypothesis test for von Mises-Fisher distribution over Kent distribution}
\alias{fishkent}
\title{
Hypothesis test for von Mises-Fisher distribution over Kent distribution
}

\description{
The null hypothesis is whether a von Mises-Fisher distribution fits the data well, where the altenrative is that Kent distribution is more suitable.
}

\usage{
fishkent(x, B = 999)
}

\arguments{
\item{x}{
A numeric matrix containing the data as unit vectors in Euclidean coordinates.
}
\item{B}{
The number of bootstrap re-samples. By default is set to 999. If it is equal to 1, no bootstrap is performed and the p-value is obtained throught the asymptotic distribution.
}
}

\details{
Essentially it is a test of rotational symmetry, whether Kent's ovalness parameter (beta) is equal to zero. This works for spherical data only.
}

\value{
A vector including:
\item{test}{
The value of the test statistic
}
\item{p-value or Bootstrap p-value}{
The p-value of the test.
}
}

\references{
Rivest L. P. (1986). Modified Kent's statistics for testing goodness of fit for the Fisher distribution in small concentrated samples. 
Statistics & Probability Letters, 4(1): 1--4.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{iagesag}, \link{vmf.mle}, \link{kent.mle}, \link{rkent}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
fishkent(x)
fishkent(x, B = 1)
iagesag(x)
}

\keyword{ Goodness of fit test }
\keyword{ von Mises-Fisher distribution }
\keyword{ Kent distribution }

