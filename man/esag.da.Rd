\name{Cross validation with ESAG discrminant analysis}
\alias{esag.da}
\title{
Cross validation for estimating the classification rate of a discrminant analysis for directional data assuming an ESAG distribution
}

\description{
Cross validation for estimating the classification rate of a discrminant analysis for directional data assuming an ESAG distribution.
}

\usage{
esag.da(y, ina, fraction = 0.2, R = 100, seed = NULL)
}

\arguments{
\item{y}{
A matrix with the data in Eulcidean coordinates, i.e. unit vectors. The matrix must have three columns, only spherical data are currently supported.
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
Tsagris M. and Alenazi A. (2019). Comparison of discriminant analysis methods on the sphere. Communications in Statistics: Case Studies, Data Analysis and Applications, 5(4), 467--491.

Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.

Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
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
esag.da(x, ina, fraction = 0.2, R = 50)
}

\keyword{ Discriminant analysis }
\keyword{ Cross validation }
\keyword{ von Mises-Fisher distribution }

