\name{Tuning of the k-NN algorithm using the arc cosinus distance}
\alias{dirknn.tune}
\title{
k-NN algorithm using the arc cosinus distance. Tuning the k neigbours
}

\description{
It estimates the percentage of correct classification via an m-fold cross valdiation. The bias is estimated as well using the algorithm suggested by Tibshirani and Tibshirani (2009) and is subtracted.
}

\usage{
dirknn.tune(z, M = 10, A = 5, ina, type = "S", mesos = TRUE, mat = NULL,
parallel = FALSE)
}

\arguments{
\item{z}{
The data, a numeric matrix with unit vectors.
}
\item{M}{
The number of folds for the m-fold cross validation, set to 10 by default.
}
\item{A}{
The maximum number of nearest neighbours, set to 5 by default. The 1 nearest neighbour is not used.
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
\item{mat}{
You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations. You can also leave it NULL and it will create folds.
}
\item{parallel}{
If you want the standard -NN algorithm to take place in parallel set this equal to TRUE.
}
}

\details{
The standard algorithm is to keep the k nearest observations and see the groups of these observations. The new observation is allocated to the most frequent seen group. The non standard algorithm is to calculate the classical mean or the harmonic mean of the k nearest observations for each group. The new observation is allocated to the group with the smallest mean distance.

We have made an eficient (not very much efficient though) memory allocation. Even if you have hundreds of thousands of observations, the computer will not clush, it will only take longer. Instead of calcualte the distance matrix once in the beginning we calcualte the distances of the out-of-sample observations from the rest. If we calculated the distance matrix in the beginning, once, the resulting matrix could have dimensions thousands by thousands. This would not fit into the memory. If you have a few hundres of observations, the runtime is about the same (maybe less, maybe more) as calculating the distance matrix in the first place.
}

\value{
A list including:
\item{per}{
The average percent of correct classification across the neighbours.
}
\item{percent}{
The bias corrected percent of correct classification.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

%\references{
%}

\author{
Michail Tsagris
R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dirknn}, \link{vmf.da}, \link{mix.vmf}
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
dirknn.tune(x, M = 5, A = 10, ina, type = "S", mesos = TRUE)
dirknn.tune(x, M = 10, A = 5, ina, type = "S", mesos = TRUE)
}

\keyword{ Directional k-NN algorithm }
\keyword{ Supervised classification }
