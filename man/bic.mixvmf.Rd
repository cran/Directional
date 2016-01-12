\name{bic.mixvmf}
\alias{bic.mixvmf}
\title{
BIC to choose the number of components in a model based clustering using mixtures of von Mises-Fisher distributions
}
\description{
BIC to choose the number of components in a model based clustering using mixtures of von Mises-Fisher distributions
}
\usage{
bic.mixvmf(x, A)
}
\arguments{
  \item{x}{
  A matrix containing directional data. 
}
  \item{A}{
  The maximum number of clusters to be tested. Default value is 3. 
}
}
\details{
If the data are not unit vectors, they are transformed into unit vectors. 
}
\value{
A list including:
\item{BIC}{
The BIC values for all the models tested. 
}
\item{A plot}{
A plot of the BIC values. 
}
}
\references{
Hornik, K. and Grun, B. (2014). movMF: An R package for fitting mixtures of von misesfisher distributions. Journal of Statistical Software, 58(10):1--31.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{mix.vmf}, \link{rmixvmf}, \link{mixvmf.contour} 
}
}
\examples{
x <- as.matrix( iris[, 1:4] )
x <- x/sqrt(rowSums(x^2))
bic.mixvmf(x, 5)  

}
\keyword{ Mixtures of von Mises-Fisher distributions}
\keyword{ BIC of the mixture models }
