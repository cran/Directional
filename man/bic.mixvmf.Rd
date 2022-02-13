\name{BIC for the model based clustering using mixtures of von Mises-Fisher distributions}
\alias{bic.mixvmf}
\title{
BIC to choose the number of components in a model based clustering using mixtures of von Mises-Fisher distributions
}

\description{
BIC to choose the number of components in a model based clustering using mixtures of von Mises-Fisher distributions
}

\usage{
bic.mixvmf(x, A, n.start = 20)
}

\arguments{
\item{x}{
A matrix containing directional data.
}
\item{A}{
The maximum number of clusters to be tested. Default value is 3.
}
\item{n.start}{
The number of random starts to try. See also R's built-in function \code{\link{kmeans}} for more information about this.
}
}

\details{
If the data are not unit vectors, they are transformed into unit vectors.
}

\value{
A plot of the BIC values and a list including:
\item{BIC}{
The BIC values for all the models tested.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
Hornik, K. and Grun, B. (2014). movMF: An R package for fitting mixtures of von Mises-Fisher distributions. Journal of Statistical Software, 58(10):1--31.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mixvmf.mle}, \link{rmixvmf}, \link{mixvmf.contour}
}
}

\examples{
x <- as.matrix( iris[, 1:4] )
x <- x / sqrt( rowSums(x^2) )
bic.mixvmf(x, 5)
}

\keyword{ Mixtures of von Mises-Fisher distributions }
\keyword{ BIC of the mixture models }
