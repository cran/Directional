\name{Prediction in discriminant analysis based on Purkayastha distribution}
\alias{purkada.pred}
\title{
Prediction of a new observation using discriminant analysis based on Purkayastha distribution
}

\description{
Prediction of a new observation using discriminant analysis based on Purkayastha distribution.
}

\usage{
purkada.pred(ynew, y, ina)
}

\arguments{
\item{ynew}{
The new observation(s) whose group is to be predicted.
A numerical vector with data expressed in radians, or a matrix with two columns (cos and sin) for circular data.
Or a matrix with 3 columns (unit vectors) for spherical data.
}
\item{y}{
A numerical vector with data expressed in radians, or a matrix with two columns (cos and sin) for circular data.
Or a matrix with 3 columns (unit vectors) for spherical data.
}
\item{ina}{
A vector indicating the groups of the data y.
}
}

\details{
Prediction of the class of a new spherical vector assuming ESAG distribution.
}

\value{
A vector with the predicted group of each new observation.
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
\code{\link{esag.da}, \link{vmfda.pred}, \link{dirknn}, \link{knn.reg}
}
}

\examples{
m1 <- rnorm(3)
m2 <- rnorm(3) + 0.5
y <- rbind( rvmf(100, m1, 3), rvmf(80, m2, 5) )
ina <- c( rep(1,100), rep(2, 80) )
ynew <- rbind(rvmf(10, m1, 10), rvmf(10, m2, 5))
id <- rep(1:2, each = 10)
g <- purkada.pred(ynew, y, ina)
table(id, g)
}

\keyword{ Discriminant analysis }
\keyword{ Von Mises-Fisher distributions }
