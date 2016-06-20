\name{vmf.da}
\alias{vmf.da}
\title{
Cross validation for estimating the classification rate of a discrminant analysis for directional data assuming a von Mises-Fisher distribution
}

\description{
Cross validation for estimating the classification rate of a discrminant analysis for directional data assuming a von Mises-Fisher distribution.
}

\usage{
vmf.da(x, ina, fraction = 0.2, R = 1000, seed = FALSE)
}

\arguments{
\item{x}{
A matrix with the dataset in Eulcidean coordinates, i.e. unit vectors.
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
If seed is TRUE, the results will always be the same.
}
}

\details{
A repeated cross validation procedure is performed to estimate the rate of correct classification.
}

\value{
A list including:
\item{percent}{
The estimated percent of correct classification and two estimated standard deviations. The one is the standard devation of the rates and the other is assuming a binomial distribution.
}
\item{ci}{
Three types of confidence intervals, the standard one, another one based on the binomial distribution and the third one is the empirical one, which calcualtes the upper and lower 2.5\% of the rates.
}
}

\references{
Morris, J. E., & Laycock, P. J. (1974). Discriminant analysis of directional data. Biometrika, 61(2): 335-341.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmfda.pred}, \link{mix.vmf}, \link{vmf}, \link{dirknn}
}
}

\examples{
x <- rvmf(100, rnorm(4), 15)
ina <- rep(1:2, each = 50)
vmf.da(x, ina, fraction = 0.2, R = 200, seed = FALSE)
}

\keyword{ Discriminant analysis }
\keyword{ Cross validation }
\keyword{ von Mises-Fisher distribution }

