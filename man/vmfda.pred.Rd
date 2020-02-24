\name{Prediction in discriminant analysis based on von Mises-Fisher distribution}
\alias{vmfda.pred}
\title{
Prediction of a new observation using discriminant analysis based on von Mises-Fisher distribution
}

\description{
Prediction of the class of a new observation using discriminant analysis based on von Mises-Fisher distribution.
}

\usage{
vmfda.pred(xnew, x, ina)
}

\arguments{
\item{xnew}{
The new observation(s) (unit vector(s)) whose group is to be predicted.
}
\item{x}{
A data matrix with unit vectors, i.e. directional data.
}
\item{ina}{
A vector indicating the groups of the data x.
}
}

\details{
Discriminant analysis assuming von Mises-Fisher distributions.
}

\value{
A vector with the predicted group of each new observation.
}

\references{
Tsagris M. and Alenazi A. (2019). Comparison of discriminant analysis methods on the sphere. Communications in Statistics: Case Studies, Data Analysis and Applications, 5(4), 467--491.

Morris J. E. and Laycock P. J. (1974). Discriminant analysis of directional data. Biometrika, 61(2): 335-341.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.da}, \link{mix.vmf}, \link{mix.vmf}, \link{dirknn}, \link{knn.reg}
}
}

\examples{
m1 <- rnorm(5)
m2 <- rnorm(5)
x <- rbind( rvmf(100, m1, 5), rvmf(80, m2, 10) )
ina <- c( rep(1,100), rep(2, 80) )
y <- rbind(rvmf(10, m1, 10), rvmf(10, m2, 5))
id <- rep(1:2, each = 10)
g <- vmfda.pred(y, x, ina)
table(id, g)
}

\keyword{ Discriminant analysis }
\keyword{ Von Mises-Fisher distributions }
