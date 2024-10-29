\name{Density of some (hyper-)spherical distributions}
\alias{dvmf}
\alias{iagd}
\alias{dpurka}
\alias{dspcauchy}
\alias{dpkbd}
\title{
Density of some (hyper-)spherical distributions
}

\description{
Density of some (hyper-)spherical distributions.
}

\usage{
dvmf(y, mu, k, logden = FALSE )
iagd(y, mu, logden = FALSE)
dpurka(y, theta, a, logden = FALSE)
dspcauchy(y, mu, rho, logden = FALSE)
dpkbd(y, mu, rho, logden = FALSE)
}

\arguments{
\item{y}{
A matrix or a vector with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{mu}{
The mean direction (unit vector) of the von Mises-Fisher, the IAG, the spherical Cauchy distribution,
or of the Poisson kernel-based distribution.
}
\item{theta}{
The mean direction (unit vector) of the Purkayastha distribution.
}
\item{k}{
The concentration parameter of the von Mises-Fisher distribution.
}
\item{a}{
The concentration parameter of the Purkayastha distribution.
}
\item{rho}{
The \eqn{\rho} parameter of the spherical Cauchy distribution, or of the
Poisson kernel-based distribution.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the von Mises-Fisher, of the IAG, of the Purkayastha, of the spherical Cauchy distribution, or of the Poisson kernel-based distribution is computed.
}

\value{
A vector with the (log) density values of y.
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through Maximum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70--83

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution. Communications in Statistics-Theory and Methods, 19(6): 1973--1986.

Kato S. and McCullagh P. (2020). Some properties of a Cauchy family on the sphere derived from the Mobius transformations. Bernoulli, 26(4): 3224--3248.
https://arxiv.org/pdf/1510.07679.pdf

Golzy M. and Markatou M. (2020). Poisson kernel-based clustering on the sphere:
convergence properties, identifiability, and a method of sampling.
Journal of Computational and Graphical Statistics, 29(4): 758--770.

Sablica L., Hornik K. and Leydold J. (2023). Efficient sampling from the PKBD distribution.
Electronic Journal of Statistics, 17(2): 2180--2209.

Zehao Yu and Xianzheng Huang (2024). A new parameterization for elliptically symmetric angular Gaussian distributions of arbitrary dimension. Electronic Journal of Statististics, 18(1): 301--334.

Tsagris M. and Papastamoulis P. (2024). Directional data analysis using the spherical Cauchy and the Poisson kernel-based distribution. https://arxiv.org/pdf/2409.03292.
}

\author{
Michail Tsagris and  Zehao Yu.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Zehao Yu \email{zehaoy@email.sc.edu}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kent.mle}, \link{rkent}, \link{esag.mle}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- rvmf(1000, m = m, k = 10)
dvmf(y, k=10, m)
}

