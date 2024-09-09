\name{Density of the spherical ESAG and Kent distributions and of the ESAG distribution in arbitrary dimensions}
\alias{desag}
\alias{dkent}
\alias{dESAGd}
\title{
Density of the spherical ESAG and Kent distributions
}

\description{
Density of the spherical ESAG and Kent distributions.
}

\usage{
desag(y, mu, gam, logden = FALSE)
dkent(y, G, param, logden = FALSE)
dESAGd(y, mu, gam, logden = FALSE)
}

\arguments{
\item{y}{
A matrix or a vector with the data expressed in Euclidean coordinates, i.e. unit vectors. For the dESAGd it can have any dimension.
}
\item{mu}{
The mean vector the ESAG distribution.
}
\item{gam}{
The two \eqn{\gamma} parameters of the ESAG distribution.
}
\item{G}{
For the Kent distribution only, a 3 x 3 matrix whose first column is the mean direction. The second and third columns are the major and minor axes respectively.
}
\item{param}{
For the Kent distribution a vector with the concentration \eqn{\kappa} and ovalness \eqn{\beta} parameters. The \eqn{\psi} has been absorbed inside the matrix G.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the spherical ESAG or Kent distribution, or of the ESAG distribution in arbitrary dimensions is computed.
}

\value{
A vector with the (log) density values of y.
}

\references{
Zehao Yu and Xianzheng Huang (2024). A new parameterization for elliptically symmetric angular Gaussian distributions of arbitrary dimension. Electronic Journal of Statististics, 18(1): 301--334.

Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.

Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society,
Series B, 44(1): 71--80.

Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris and Zehao Yu.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Zehao Yu \email{Zzehaoy@email.sc.edu}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kent.mle}, \link{rkent}, \link{esag.mle}
}
}

\examples{
m <- colMeans( as.matrix( iris[, 1:3] ) )
y <- rkent(1000, k = 10, m = m, b = 4)
mod <- kent.mle(y)
dkent( y, G = mod$G, param = mod$param )
}

