\name{The k-nearest neighbours using the cosinus distance}
\alias{cosnn}
\title{
The k-nearest neighbours using the cosinus distance
}

\description{
The k-nearest neighbours using the cosinus distance.
}

\usage{
cosnn(xnew, x, k = 5, index = FALSE, rann = FALSE)
}

\arguments{
\item{xnew}{
The new data whose k-nearest neighbours are to be found.
}
\item{x}{
The data, a numeric matrix with unit vectors.
}
\item{k}{
The number of nearest neighbours, set to 5 by default. It can also be a vector with many values.
}
\item{index}{
If you want the indices of the closest observations set this equal to TRUE.
}
\item{rann}{
If you have large scale datasets and want a faster k-NN search, you can use kd-trees implemented in the R package "RANN". In this case you must set this argument equal to TRUE.
}
}

\details{
The shortest distances or the indices of the k-nearest neighbours using the cosinus distance are returned.}

\value{
A matrix with the shortest distance of each xnew from x, if index is FALSE, or the indices of the nearest neighbours of each xnew from x, if index is TRUE.
}

\references{
Tsagris M. and Alenazi A. (2019). Comparison of discriminant analysis methods on the sphere. Communications in Statistics: Case Studies, Data Analysis and Applications, 5(4), 467--491.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dirknn}, \link{dirknn.tune}
}
}

\examples{
xnew <- rvmf(10, rnorm(3), 5)
x <- rvmf(50, rnorm(3), 5)
a <- cosnn(xnew, x, k = 5)
b <- cosnn(xnew, x, k = 5, index = TRUE)
}
