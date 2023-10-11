\name{Simulation of random values from rotationally symmetric distributions}
\alias{rvmf}
\alias{riag}
\alias{rspcauchy}
\alias{rpkbd}
\title{
Simulation of random values from rotationally symmetric distributions
}

\description{
Simulation of random values from rotationally symmetric distributions. The data can be spherical or hyper-spherical.
}

\usage{
rvmf(n, mu, k)
riag(n, mu)
rspcauchy(n, mu, rho)
rpkbd(n, mu, rho)
}

\arguments{
\item{n}{
The sample size.
}
\item{mu}{
A unit vector showing the mean direction for the von Mises-Fisher or the spherical Cauchy distribution. The mean vector of the Independent Angular Gaussian distribution does not have to be a unit vector.
}
\item{k}{
The concentration parameter (\eqn{\kappa}) of the von Mises-Fisher distribution. If \eqn{\kappa=0}, random values from the spherical uniform will be drwan.
}
\item{rho}{
The \eqn{\rho} parameter of the spherical Cauchy or the Poisson kernel based distribution.
}
}

\details{
The von Mises-Fisher uses the rejection smapling suggested by Wood (1994). For the Independent Angular Gaussian, values are generated from a multivariate normal distribution with the given mean vector and the identity matrix as the covariance matrix. Then each vector becomes a unit vector. For the spherical Cauchy distribution the algortihm is described in Kato and McCullagh (2020) and for the Poisson kernel based distribution, it is described in Sablica, Hornik and Leydold (2023).
}

\value{
A matrix with the simulated data.
}

\references{
Wood A.T.A. (1994). Simulation of the von Mises Fisher distribution.
Communications in Statistics-Simulation and Computation, 23(1): 157--164.

Dhillon I. S. and Sra S. (2003). Modeling data using directional distributions.
Technical Report TR-03-06, Department of Computer Sciences, The University of Texas at Austin.
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.75.4122&rep=rep1&type=pdf

Kato S. and McCullagh P. (2020). Some properties of a Cauchy family on the sphere derived from the Mobius transformations. Bernoulli, 26(4): 3224--3248.
https://arxiv.org/pdf/1510.07679.pdf

Sablica L., Hornik K. and Leydold J. (2023). Efficient sampling from the PKBD distribution.
Electronic Journal of Statistics, 17(2): 2180--2209.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.mle}, \link{iag.mle} \link{rfb}, \link{racg}, \link{rvonmises}, \link{rmixvmf}
}
}

\examples{
m <- rnorm(4)
m <- m/sqrt(sum(m^2))
x <- rvmf(100, m, 25)
m
vmf.mle(x)
}

\keyword{ Random values simulation }
\keyword{ von Mises-Fisher distribution }
