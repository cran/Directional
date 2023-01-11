\name{Rayleigh's test of uniformity}
\alias{rayleigh}
\title{
Rayleigh's test of uniformity
}

\description{
It checkes whether the data are uniformly distributed on the sphere or hypersphere.
}

\usage{
rayleigh(x, modif = TRUE, B = 999)
}

\arguments{
\item{x}{
A matrix containing the data, unit vectors.
}
\item{modif}{
If modif is TRUE, the modification as suggested by Jupp (2001) is used.
}
\item{B}{
If B is greater than 1, bootstap calibation os performed. If it is equal to 1, classical theory is used.
}
}

\details{
The Rayleigh test of uniformity is not the best, when there are two antipodal mean directions. In this case it will fail. It is good to test whether there is one mean direction or not. To put it differently, it tests whether the concentration parameter of the Fisher distribution is zero or not.
}

\value{
A vector including:
\item{test}{
The value of the test statistic.
}
\item{p-value or Bootstrap p-value}{
The (bootstrap) p-value of the test.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Jupp, P. E. (2001). Modifications of the rayleigh and bingham tests for uniformity of directions. Journal of Multivariate Analysis, 77(2):1-20.

Rayleigh, L. (1919). On the problem of random vibrations, and of random flights in one, two, or three dimensions. The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science, 37(220):321-347.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.mle}, \link{meandir.test}, \link{acg.mle}
}
}

\examples{
x <- rvmf(100, rnorm(5), 1)  ## Fisher distribution with low concentration
rayleigh(x)
}

\keyword{ Uniform distribution on the (hyper)-sphere }
\keyword{ Hypothesis test }
