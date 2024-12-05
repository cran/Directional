\name{BIC for the model based clustering using mixtures of rotationally symmetric distributions}
\alias{bic.mixvmf}
\alias{bic.mixspcauchy}
\alias{bic.mixpkbd}
\title{
BIC to choose the number of components in a model based clustering using mixtures of rotationally symmetric distributions
}

\description{
BIC to choose the number of components in a model based clustering using mixtures of rotationally symmetric distributions
}

\usage{
bic.mixvmf(x, G = 5, n.start = , tol = 1e-6, maxiters = 500)
bic.mixspcauchy(x, G = 5, n.start = 5, tol = 1e-6, maxiters = 500)
bic.mixpkbd(x, G = 5, n.start = 5, tol = 1e-6, maxiters = 500)
}

\arguments{
\item{x}{
A matrix containing directional data.
}
\item{G}{
The maximum number of clusters to be tested. Default value is 5.
}
\item{n.start}{
The number of random starts to try. See also R's built-in function \code{\link{kmeans}} for more information about this.
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxiters}{
The maximum number of iterations to perform.
}
}

\details{
The function computes the BIC (and ICL) to decide on the optimal number of clusters when using mixtures of von Mises-Fisher, mixtures of spherical Cauchy or mixtures of Poisson kernel-based distributions.
}

\value{
A plot of the BIC values and a list including:
\item{bic}{
The BIC values for all the models tested.
}
\item{icl}{
The ICL values for all the models tested.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
Hornik, K. and Grun, B. (2014). movMF: An R package for fitting mixtures of von Mises-Fisher distributions. Journal of Statistical Software, 58(10): 1--31.

Biernacki C., Celeux G. and Govaert, G. (2000). Assessing a mixture model for clustering with the integrated completed likelihood. IEEE Transactions on Pattern Analysis and Machine Intelligence, 22(7): 719--725.

Tsagris M., Papastamoulis P. and Kato S. (2024). Directional data analysis using the spherical Cauchy and the Poisson kernel-based distribution. https://arxiv.org/pdf/2409.03292
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
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
bic.mixvmf(x)
}

