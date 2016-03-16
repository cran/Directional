\name{knn.reg}
\alias{knn.reg}

\title{
k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables
}

\description{
k-NN regression with Euclidean or (hyper-)spherical response and or predictor variables.
}

\usage{
knn.reg(xnew, y, x, k = 5, res = "eucl", type = "euclidean", estim = "arithmetic")
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
The number of nearest neighbours, set to 5 by default.
}
\item{res}{
The type of the response variable. If it is Euclidean, set this argument equal to "res".  If it is a unit vector set it to res="spher".
}
\item{type}{
The type of distance to be used. This determines the nature of the predictor variables. This is actually the argument "method" of the distance function in R. The default value is "euclidean". If you use the Euclidean distance, the package "fields" is used. The "rdist" function of that packages is about 3 times faster than the standard built-in "dist". R has several options the type of the distance. Just type ?dist in R and see the methods. Any method can be given here. If you have unit vectors in general, you should put type="spher", so that the cosinus distance is calculated.
}
\item{estim}{
Once the k observations whith the smallest distance are discovered, what should the prediction be? The arithmetic average of the corresponding y values be used estim="arithmetic" or their harmonic average estim="harmonic".
}
}

\details{
This funciton covers a broad range of data, Euclidean and spherical, along with their combinations.
}

\value{
A matrix (or a vector in the case of Euclidean data) with the predicted response values.
}


%\references{
%
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
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
x <- iris[, 2:4]
x <- x/ sqrt( rowSums(x^2) )  ## Euclidean response and spherical predictors
a <- knn.reg(x, y, x, k = 5, res = "eucl", type = "spher", estim = "arithmetic")

y <- iris[, 2:4]
y <- y/ sqrt( rowSums(y^2) )  ## Spherical response and Euclidean predictor
x <- iris[, 1]
a <- knn.reg(x, y, x, k = 5, res = "spher", type = "euclidean", estim = "arithmetic")
}

\keyword{ k-NN regression }
\keyword{ Euclidean data }
\keyword{ Spherical data }
