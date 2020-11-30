\name{Tuning of the k-NN regression}
\alias{knnreg.tune}
\title{
Tuning of the k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables
}

\description{
Tuning of the k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables. It estimates the percentage of correct classification via an m-fold cross valdiation. The bias is estimated as well using the algorithm suggested by Tibshirani and Tibshirani (2009) and is subtracted.
}

\usage{
knnreg.tune(y, x, nfolds = 10, A = 10, ncores = 1, res = "eucl",
estim = "arithmetic", folds = NULL, seed = FALSE, graph = FALSE)
}

\arguments{
\item{y}{
The currently available data, the response variables values. A matrix with either euclidean (univariate or multivariate) or (hyper-)spherical data. If you have a circular response, say u, transform it to a unit vector via (cos(u), sin(u)).
}
\item{x}{
The currently available data, the predictor variables values. A matrix with either euclidean (univariate or multivariate) or (hyper-)spherical data. If you have a circular response, say u, transform it to a unit vector via (cos(u), sin(u)).
}
\item{nfolds}{
How many folds to create?
}
\item{A}{
The maximum number of nearest neighbours, set to 5 by default, starting from the 1 nearest neighbor.
}
\item{ncores}{
How many cores to use. This is taken into account only when the predictor variables are spherical.
}
\item{res}{
The type of the response variable. If it is Euclidean, set this argument equal to "res".  If it is a unit vector set it to res="spher".
}
\item{estim}{
Once the k observations whith the smallest distance are discovered, what should the prediction be? The arithmetic average of the corresponding y values be used estim="arithmetic" or their harmonic average estim="harmonic".
}
\item{folds}{
Do you already have a list with the folds? If not, leave this NULL.
}
\item{seed}{
If seed is TRUE, the results will always be the same.
}
\item{graph}{
If this is TRUE a graph with the results will appear.
}
}

\details{
Tuning of the k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables. It estimates the percentage of correct classification via an m-fold cross valdiation. The bias is estimated as well using the algorithm suggested by Tibshirani and Tibshirani (2009) and is subtracted. The sum of squares of prediction is used in the case of Euclidean responses. In the case of spherical responses the \eqn{\sum_{\hat{y}_i^T}y_i} is calculated.
}

\value{
A list including:
\item{crit}{
The value of the criterion to minimise/maximise for all values of the nearest neighbours.
}
\item{best_k}{
The best value of the nearest neighbours.
}
\item{performance}{
The bias corrected optimal value of the criterion, along wit the estimated bias. For the case of Euclidean reponse this will be higher than the crit and for the case or spherical responses it will be lower than crit.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

%\references{
%}

\author{
Michail Tsagris
R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{knn.reg}, \link{spher.reg}, \link{dirknn.tune}
}
}

\examples{
y <- iris[, 1]
x <- iris[, 2:4]
x <- x/ sqrt( rowSums(x^2) )  ## Euclidean response and spherical predictors
knnreg.tune(y, x, A = 5, res = "eucl", estim = "arithmetic")

y <- iris[, 1:3]
y <- y/ sqrt( rowSums(y^2) )  ## Spherical response and Euclidean predictor
x <- iris[, 2]
knnreg.tune(y, x, A = 5, res = "spher", estim = "arithmetic")
}


\keyword{ k-NN regression }
\keyword{ Cross-validation }
\keyword{ Euclidean data }
\keyword{ Spherical data }
