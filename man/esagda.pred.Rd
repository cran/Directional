\name{Prediction in discriminant analysis based on ESAG distribution}
\alias{esagda.pred}
\title{
Prediction of a new observation using discriminant analysis based on ESAG distribution
}

\description{
Prediction of a new observation using discriminant analysis based on ESAG distribution.
}

\usage{
esagda.pred(ynew, y, ina)
}

\arguments{
\item{ynew}{
The new observation(s) (unit vector(s)) whose group is to be predicted.
}
\item{y}{
A data matrix with unit vectors, i.e. spherical directional data.
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
Tsagris M. and Alenazi A. (2019). Comparison of discriminant analysis methods on the sphere. Communications in Statistics: Case Studies, Data Analysis and Applications, 5(4), 467--491.

Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2017). An Elliptically Symmetric Angular
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
g <- esagda.pred(ynew, y, ina)
table(id, g)
}

\keyword{ Discriminant analysis }
\keyword{ Von Mises-Fisher distributions }
