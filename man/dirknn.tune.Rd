\name{Tuning of the k-NN algorithm using the arc cosinus distance}
\alias{dirknn.tune}
\title{
k-NN algorithm using the arc cosinus distance. Tuning the k neigbours
}

\description{
It estimates the percentage of correct classification via an m-fold cross valdiation. The bias is estimated as well using the algorithm suggested by Tibshirani and Tibshirani (2009) and is subtracted.
}

\usage{
dirknn.tune(z, nfolds = 10, A = 5, ina, type = "S", mesos = TRUE, folds = NULL,
parallel = FALSE, stratified = TRUE, seed = FALSE)
}

\arguments{
\item{z}{
The data, a numeric matrix with unit vectors.
}
\item{nfolds}{
How many folds to create?
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
\item{folds}{
Do you already have a list with the folds? If not, leave this NULL.
}
\item{parallel}{
If you want the standard -NN algorithm to take place in parallel set this equal to TRUE.
}
\item{stratified}{
Should the folds be created in a stratified way? i.e. keeping the distribution of the groups similar through all folds?
}
\item{seed}{
If seed is TRUE, the results will always be the same.
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
R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
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
dirknn.tune(x, nfolds = 5, A = 10, ina, type = "S", mesos = TRUE)
dirknn.tune(x, nfolds = 10, A = 5, ina, type = "S", mesos = TRUE)
}

\keyword{ Directional k-NN algorithm }
\keyword{ Supervised classification }
