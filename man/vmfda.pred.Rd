\name{vmfda.pred}
\alias{vmfda.pred}
\title{
Prediction of a new observation using discriminant analysis based on von Mises-Fisher distributions.
}
\description{
Prediction of a new observation using discriminant analysis based on von Mises-Fisher distributions.
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
Morris, J. E., & Laycock, P. J. (1974). Discriminant analysis of directional data. Biometrika, 61(2), 335-341.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmf.da}, \link{mix.vmf}}
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
