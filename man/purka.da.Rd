\name{Cross validation with Purkayastha discrminant analysis}
\alias{purka.da}
\title{
Cross validation for estimating the classification rate of a discrminant analysis for directional data assuming a Purkayastha distribution
}

\description{
Cross validation for estimating the classification rate of a discrminant analysis for directional data assuming a Purkayastha distribution.
}

\usage{
purka.da(y, ina, fraction = 0.2, R = 100, seed = NULL)
}

\arguments{
\item{y}{
A numerical vector with data expressed in radians, or a matrix with two columns (cos and sin) for circular data.
Or a matrix with 3 columns (unit vectors) for spherical data.
}
\item{ina}{
A variable indicating the groupings.
}
\item{fraction}{
The fraction of data to be used as test set.
}
\item{R}{
The number of repetitions.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
}

\details{
A repeated cross validation procedure is performed to estimate the rate of correct classification.
}

\value{
A list including:
\item{percent}{
The estimated percent of correct classification and two estimated standard deviations. The one is the standard devation of the rates and the other is
assuming a binomial distribution.
}
\item{ci}{
Three types of confidence intervals, the standard one, another one based on the binomial distribution and the third one is the empirical one,
which calcualtes the upper and lower 2.5\% of the rates.
}
}

\references{
Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through Maximum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70-83

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution. Communications in Statistics-Theory and Methods, 19(6): 1973-1986.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.da}, \link{vmfda.pred}, \link{dirknn}
}
}

\examples{
x <- rvmf(100, rnorm(3), 15)
ina <- rep(1:2, each = 50)
purka.da(x, ina, fraction = 0.2, R = 50)
}

\keyword{ Discriminant analysis }
\keyword{ Cross validation }
\keyword{ von Mises-Fisher distribution }

