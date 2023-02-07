\name{Contour plot (on the sphere) of a mixture of von Mises-Fisher distributions}
\alias{spher.mixvmf.contour}
\title{
Contour plot (on the sphere) of a mixture of von Mises-Fisher distributions
}

\description{
The contour plot (on the sphere) of a mixture of von Mises-Fisher distributions is produced.
}

\usage{
spher.mixvmf.contour(probs, mu, k, bgcol = "snow", dat = NULL, col = NULL,
lat = 50, long = 50)
}

\arguments{
\item{probs}{
This is avector with the mixing probability of each group.
}
\item{mu}{
A matrix with the mean direction of each group.
}
\item{k}{
A vector with the concentration parameter of each group.
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
The goal of this function is for the user to see how the mixtures of von Mises-Fisher look like.
}

\value{
A plot containing the contours of the mixture distribution.
}

\references{
Kurt Hornik and  Bettina Grun (2014). movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf

Mardia K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra S. (2012). A short note on parameter approximation for von Mises-Fisher distributions:
and a fast implementation of \eqn{I_s(x)}. Computational Statistics, 27(1): 177--190.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spher.esag.contour}, \link{spher.vmf.contour}, \link{mixvmf.mle}
}
}

\examples{
\donttest{
k <- runif(3, 4, 20)
probs <- c(0.2, 0.5, 0.3)
mu <- matrix(rnorm(9, 0, 0.5), ncol = 3)
mu <- mu / sqrt( rowSums(mu^2) )
## the lat and long are decreased to 10. Increase them back to 50 to
## see the difference
spher.mixvmf.contour(probs, mu, k, lat = 10, long = 10)
}
}

