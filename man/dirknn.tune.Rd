\name{dirknn.tune}
\alias{dirknn.tune}
\title{
k-NN algorithm using the arc cosinus distance. Tuning the k neigbours.
}
\description{
It estimates the percentage of correct classification via an m-fold cross valdiation. The bias is estimated as well using the algorithm suggested by Tibshirani and Tibshirani (2009) and is subtracted. 
}
\usage{
dirknn.tune(x, M = 10, A = 5, ina, type = "S", mesos = TRUE, seed = TRUE) 
}
\arguments{
  \item{x}{
  The data, a numeric matrix with unit vectors.
}
  \item{M}{
  The number of folds for the m-fold cross validation, set to 10 by default. 
}
  \item{A}{
  The maximum number of nearest neighbours, set to 5 by default. The maximum is actually A + 1, since 1 nearest neighbour is not used. 
}
  \item{ina}{
  A variable indicating the groups of the data x.
}
  \item{type}{
  If type is "S", the standard k-NN algorithm is to be used, else "NS" for the non standard one. See below (details) for more information.
}
  \item{mesos}{
  A boolean variable used only in the case of the non standard algorithm (type="NS"). Should the average of the distances be calculated (TRUE) or not (FALSE)? If it is FALSE, the harmonic mean is calculated.
}
  \item{seed}{
  A boolean variable. If seed is TRUE then the same folds will always be produced.
}
}
\details{
The standard algorithm is to keep the k nearest observations and see the groups of these observations. The new observation is allocated to the most frequent seen group. The non standard algorithm is to calculate the classical mean or the harmonic mean of the k nearest observations for each group. The new observation is allocated to the group with the smallest mean distance. The estimated bias is calculated as Tibshirani and Tibshirani (2009) suggested.  
}
\value{
A vector including:
\item{g}{
The predicted group.
}
}
\references{
course webpage of Howard E. Haber. 
http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf

Tibshirani, Ryan J., and Robert Tibshirani. "A bias correction for the minimum error rate in cross-validation." The Annals of Applied Statistics (2009), 3(2): 822-829.
}
\author{
Michail Tsagris 
R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{dirknn}, \link{vmf.da}
}
}
\examples{
k <- runif(4, 4, 20)
prob <- c(0.2, 0.4, 0.3, 0.1)
mu <- matrix(rnorm(16), ncol = 4)
mu <- mu / sqrt( rowSums(mu^2) )
da <- rmixvmf(200, prob, mu, k) 
x <- da$x
ina <- da$id
dirknn.tune(x, M = 5, A = 10, ina, type = "S", mesos = TRUE, seed = TRUE) 
dirknn.tune(x, M = 10, A = 5, ina, type = "S", mesos = TRUE, seed = TRUE) 
}
\keyword{ Directional k-NN algorithm }
\keyword{ Supervised classification }
