\name{k-NN algorithm using the arc cosinus distance}
\alias{dirknn}
\title{
k-NN algorithm using the arc cosinus distance
}

\description{
It classifies new observations to some known groups via the k-NN algorithm.
}

\usage{
dirknn(xnew, ina, x, k = 5, mesos = TRUE, parallel = FALSE, rann = FALSE)
}

\arguments{
\item{xnew}{
The new data whose membership is to be predicted, a numeric matrix with unit vectors.
}
\item{ina}{
A variable indicating the groups of the data x.
}
\item{x}{
The data, a numeric matrix with unit vectors.
}
\item{k}{
The number of nearest neighbours, set to 5 by default. It can also be a vector with many values.
}
\item{mesos}{
A boolean variable used only in the case of the non standard algorithm (type="NS"). Should the average of the distances be calculated (TRUE) or not (FALSE)? If it is FALSE, the harmonic mean is calculated.
}
\item{parallel}{
If you want the standard -NN algorithm to take place in parallel set this equal to TRUE.
}
\item{rann}{
If you have large scale datasets and want a faster k-NN search, you can use kd-trees implemented in the R package "RANN". In this case you must set this argument equal to TRUE.
}
}

\details{
The standard algorithm is to keep the k nearest observations and see the groups of these observations. The new observation is allocated to the most frequent seen group. The non standard algorithm is to calculate the classical mean or the harmonic mean of the k nearest observations for each group. The new observation is allocated to the group with the smallest mean distance.
}
\value{
A vector including:
\item{g}{
A matrix with the predicted group(s). It has as many columns as the values of k.
}
}

\references{
Tsagris M. and Alenazi A. (2019). Comparison of discriminant analysis methods on the sphere. Communications in Statistics: Case Studies, Data Analysis and Applications, 5(4), 467--491.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dirknn.tune}, \link{vmfda.pred}, \link{mixvmf.mle}
}
}

\examples{
k <- runif(4, 4, 20)
prob <- c(0.2, 0.4, 0.3, 0.1)
mu <- matrix(rnorm(16), ncol = 4)
mu <- mu / sqrt( rowSums(mu^2) )
da <- rmixvmf(200, prob, mu, k)
nu <- sample(1:200, 180)
x <- da$x[nu, ]
ina <- da$id[nu]
xx <- da$x[-nu, ]
id <- da$id[-nu]
a1 <- dirknn(xx, ina, x, k = 5, mesos = TRUE)
a2 <- dirknn(xx, ina, x, k = 5, mesos = FALSE)
b <- vmfda.pred(xx, x, ina)
table(id, a1)
table(id, a2)
}

\keyword{ Directional k-NN algorithm }
\keyword{ Supervised classification }
