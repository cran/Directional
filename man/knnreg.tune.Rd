\name{Tuning of the k-NN regression}
\alias{knnreg.tune}
\title{
Tuning of the k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables
}

\description{
Tuning of the k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables. It estimates the percentage of correct classification via an m-fold cross valdiation. The bias is estimated as well using the algorithm suggested by Tibshirani and Tibshirani (2009) and is subtracted.
}

\usage{
knnreg.tune(y, x, M = 10, A = 10, ncores = 1, res = "eucl", type = "euclidean",
estim = "arithmetic", mat = NULL, graph = FALSE)
}

\arguments{
\item{y}{
The currently available data, the response variables values. A matrix with either euclidean (univariate or multivariate) or (hyper-)spherical data. If you have a circular response, say u, transform it to a unit vector via (cos(u), sin(u)).
}
\item{x}{
The currently available data, the predictor variables values. A matrix with either euclidean (univariate or multivariate) or (hyper-)spherical data. If you have a circular response, say u, transform it to a unit vector via (cos(u), sin(u)).
}
  \item{M}{
  The number of folds for the m-fold cross validation, set to 10 by default.
}
\item{A}{
The maximum number of nearest neighbours, set to 5 by default. The maximum is actually A + 1, since 1 nearest neighbour is not used.
}
\item{ncores}{
How many cores to use. This is taken into account only when the predictor variables are spherical.
}
\item{res}{
The type of the response variable. If it is Euclidean, set this argument equal to "res".  If it is a unit vector set it to res="spher".
}
\item{type}{
The type of distance to be used. This determines the nature of the predictor variables. This is actually the argument "method" of the distance function in R. The default value is "euclidean". R has several options the type of the distance. Just type ?dist in R and see the methods. Any method can be given here. If you have unit vectors in general, 
you should put type="spher", so that the cosinus distance is calculated.
}
\item{estim}{
Once the k observations whith the smallest distance are discovered, what should the prediction be? The arithmetic average of the corresponding y values be used estim="arithmetic" or their harmonic average estim="harmonic".
}
\item{mat}{
You can specify your own folds by giving a mat, where each column is a fold. Each column contains indices of the observations.     You can also leave it NULL and it will create folds.
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

\references{
Tibshirani, Ryan J., and Robert Tibshirani. A bias correction for the minimum error rate in cross-validation. The Annals of Applied Statistics (2009), 3(2): 822-829.
}

\author{
Michail Tsagris
R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
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
knnreg.tune(y, x, A = 5, res = "eucl", type = "spher", estim = "arithmetic",
mat = NULL, graph = TRUE)

y <- iris[, 1:3]
y <- y/ sqrt( rowSums(y^2) )  ## Spherical response and Euclidean predictor
x <- iris[, 2]
knnreg.tune(y, x, A = 5, res = "spher", type = "euclidean", estim = "arithmetic",
mat = NULL, graph = TRUE)
}


\keyword{ k-NN regression }
\keyword{ Cross-validation }
\keyword{ Euclidean data }
\keyword{ Spherical data }
