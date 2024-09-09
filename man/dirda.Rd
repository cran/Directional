\name{Prediction in discriminant analysis based on some distributions}
\alias{dirda}
\title{
Prediction of a new observation using discriminant analysis based on some distributions
}

\description{
Prediction of a new observation using discriminant analysis based on some distributions.
}

\usage{
dirda(xnew, x, ina, type = c("vmf", "iag", "esag", "kent", "sc", "pkbd", "purka") )
}

\arguments{
\item{xnew}{
The new observation(s) (unit vector(s)) whose group is to be predicted.
}
\item{x}{
A data matrix with unit vectors, i.e. spherical directional data.
}
\item{ina}{
A vector indicating the groups of the data y.
}
\item{type}{
The type of classifier to use. The avaliable options are "vmf" (von Mises-Fisher distribution),
"iag" (IAG distribution), "esag" (ESAG distribution), "kent" (Kent distribution), "sc" and "sc2" (spherical Cauchy distribution), "pkbd" and "pkbd2" (Poisson kernel based distribution), and "purka" (Purkayastha distribution). The difference between "sc" and "sc2" and between "pkbd" and "pkbd2" is that the first uses the Newton-Raphson algorithm and it is faster, whereas the second uses a hybrid algorithm that does not require the Hessian matrix, but in large dimensions the second will be faster. You can chose any of them or all of them. Note that "kent" works only with spherical data.
}
}

\details{
Prediction of the class of a new (hyper-)spherical vector assuming some distributions.
}

\value{
A vector with the predicted group of each new observation.
}

\references{
Tsagris M. (2024). Directional data analysis using the spherical Cauchy and the Poisson-kernel based distribution.
https://arxiv.org/pdf/2409.03292

Tsagris M. and Alenazi A. (2019). Comparison of discriminant analysis methods on the sphere.
Communications in Statistics: Case Studies, Data Analysis and Applications, 5(4): 467--491.

Morris J. E. and Laycock P. J. (1974). Discriminant analysis of directional data. Biometrika, 61(2): 335--341.

Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{dirda.cv}, \link{vm.nb}, \link{dirknn}, \link{knn.reg}
}
}

\examples{
m1 <- rnorm(3)
m2 <- rnorm(3) + 0.5
x <- rbind( rvmf(100, m1, 3), rvmf(80, m2, 5) )
ina <- c( rep(1,100), rep(2, 80) )
xnew <- rbind(rvmf(10, m1, 10), rvmf(10, m2, 5))
id <- rep(1:2, each = 10)
g <- dirda(xnew, x, ina, type = "vmf")
table(id, g[, 1])
}

\keyword{ Discriminant analysis }
\keyword{ Von Mises-Fisher distributions }
