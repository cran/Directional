\name{k-NN regression}
\alias{knn.reg}

\title{
k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables
}

\description{
k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables.
}

\usage{
knn.reg(xnew, y, x, k = 5, res = "eucl", estim = "arithmetic")
}

\arguments{
\item{xnew}{
The new data, new predictor variables values. A matrix with either euclidean (univariate or multivariate) or (hyper-)spherical data. If you have a circular response, say u, transform it to a unit vector via (cos(u), sin(u)). If xnew = x, you will get the fitted values.
}
\item{y}{
The currently available data, the response variables values. A matrix with either euclidean (univariate or multivariate) or (hyper-)spherical data. If you have a circular response, say u, transform it to a unit vector via (cos(u), sin(u)).
}
\item{x}{
The currently available data, the predictor variables values. A matrix with either euclidean (univariate or multivariate) or (hyper-)spherical data. If you have a circular response, say u, transform it to a unit vector via (cos(u), sin(u)).
}
\item{k}{
The number of nearest neighbours, set to 5 by default. This can also be a vector with many values.
}
\item{res}{
The type of the response variable. If it is Euclidean, set this argument equal to "res".  If it is a unit vector set it to res="spher".
}
\item{estim}{
Once the k observations whith the smallest distance are discovered, what should the prediction be? The arithmetic average of the corresponding y values be used estim="arithmetic" or their harmonic average estim="harmonic".
}
}

\details{
This function covers a broad range of data, Euclidean and spherical, along with their combinations.
}

\value{
A list with as many elements as the number of values of k. Each element in the list contains a matrix (or a vector in the case of Euclidean data) with the predicted response values.
}


%\references{
%
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{knnreg.tune}, \link{spher.reg}, \link{spml.reg}
}
}

\examples{
y <- iris[, 1]
x <- as.matrix(iris[, 2:4])
x <- x/ sqrt( rowSums(x^2) )  ## Euclidean response
a <- knn.reg(x, y, x, k = 5, res = "eucl", estim = "arithmetic")

y <- iris[, 2:4]
y <- y / sqrt( rowSums(y^2) ) ## Spherical response
x <- iris[, 1]
a <- knn.reg(x, y, x, k = 5, res = "spher", estim = "arithmetic")
}

\keyword{ k-NN regression }
\keyword{ Euclidean data }
\keyword{ Spherical data }
