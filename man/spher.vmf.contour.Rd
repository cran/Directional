\name{Contour plot (on the sphere) of some spherical rotationally symmetric distributions}
\alias{spher.vmf.contour}
\alias{spher.purka.contour}
\alias{spher.spcauchy.contour}
\alias{spher.pkbd.contour}
\title{
Contour plot (on the sphere) of some spherical rotationally symmetric distributions
}

\description{
The contour plot (on the sphere) of some spherical rotationally symmetric distributions is produced.
}

\usage{
spher.vmf.contour(mu, k, bgcol = "snow", dat = NULL, col = NULL,
lat = 50, long = 50)
spher.purka.contour(theta, a, bgcol = "snow", dat = NULL, col = NULL,
lat = 50, long = 50)
spher.spcauchy.contour(mu, rho, bgcol = "snow", dat = NULL, col = NULL,
lat = 50, long = 50)
spher.pkbd.contour(mu, rho, bgcol = "snow", dat = NULL, col = NULL,
lat = 50, long = 50)
}

\arguments{
\item{mu}{
The mean or the median direction, depending on the distribution, a unit vector.
}
\item{theta}{
The mean direction (unit vector) of the Purkayastha distribution.
}
\item{k}{
The concentration parameter (\eqn{\kappa}) of the von Mises-Fisher distribution.
}
\item{a}{
The concentration parameter (\eqn{\alpha}) of the Purkayastha distribution.
}
\item{rho}{
The concentration parameter (\eqn{\rho}) of the spherical Cauchy distribution.
}
\item{bgcol}{
The color of the surface of the sphere.
}
\item{dat}{
If you have you want to plot supply them here. This has to be a numerical matrix with three columns, i.e. unit vectors.
}
\item{col}{
If you supplied data then choose the color of the points. If you did not choose a color, the points will appear in red.
}
\item{lat}{
A positive number determing the range of degrees to move left and right from the latitude center. See the example to better understand this argument.
}
\item{long}{
A positive number determing the range of degrees to move up and down from the longitude center. See the example to better understand this argument.
}
}

\details{
The goal of this function is for the user to see how the von Mises-Fisher, the Purkayastha, the spherical Cauchy or the Poisson kernel-based distribution looks like.
}

\value{
A plot containing the contours of the distribution.
}

\references{
Tsagris M., Papastamoulis P. and Kato S. (2024). Directional data analysis using the spherical Cauchy and the Poisson kernel-based distribution. https://arxiv.org/pdf/2409.03292.

Mardia K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra S. (2012). A short note on parameter approximation for von Mises-Fisher distributions:
and a fast implementation of \eqn{I_s(x)}. Computational Statistics, 27(1): 177--190.

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through
Maximum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70--83.

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution.
Communications in Statistics-Theory and Methods, 19(6): 1973--1986.

Kato S. and McCullagh P. (2020). Some properties of a Cauchy family on the sphere derived from the Mobius transformations. Bernoulli, 26(4): 3224--3248. https://arxiv.org/pdf/1510.07679.pdf

Golzy M. and Markatou M. (2020). Poisson kernel-based clustering on the sphere:
convergence properties, identifiability, and a method of sampling.
Journal of Computational and Graphical Statistics, 29(4): 758--770.

Sablica L., Hornik K. and Leydold J. (2023). Efficient sampling from the PKBD distribution.
Electronic Journal of Statistics, 17(2): 2180--2209.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spher.esag.contour}, \link{spher.mixvmf.contour}, \link{kent.contour}
}
}

\examples{
\donttest{
mu <- colMeans( as.matrix( iris[, 1:3] ) )
mu <- mu / sqrt( sum(mu^2) )
## the lat and long are decreased to 30. Increase them back to 50 to
## see the difference
spher.spcauchy.contour(mu, 0.7, lat = 30, long = 30)
}
}

