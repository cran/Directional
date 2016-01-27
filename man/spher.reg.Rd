\name{spher.reg}
\alias{spher.reg}
\title{
Spherical-Spherical regression.
}
\description{
It performs regression when both the dependent and independent variables are spherical.
}
\usage{
spher.reg(y, x, euclidean = FALSE, rads = FALSE)
}
\arguments{
  \item{y}{
  The dependent variable; a matrix with either two columns, latitude and longitude, either in radians or in degrees. Alternatively it is a matrix with three columns, unit vectors.
}
  \item{x}{
    The dependent variable; a matrix with either two columns, latitude and longitude, either in radians or in degrees. Alternatively it is a matrix with three columns, unit vectors. The two matrices must agree in the scale and dimensions.
}
  \item{euclidean}{
Are the data in Eulcidean coordinates, are they unit vectors. If yes, then this argument is TRUE and FALSE otherwise.
}
  \item{rads}{
If the data are expressed in latitude and longitude then it matter to know if they are in radians or degrees. If they are in radians, then this should be TRUE and FALSE otherwise. If the previous argument, euclidean, is TRUE, this one does not matter what its value is.
}
}
\details{
Spherical regression as proposed by Chang (1986) is implemented. If the estimated rotation matrix has a determinant equal to -1, singualr value decomposition is performed and the last unit vector of the second matrix is multiplied by -1.
}
\value{
A list including:
\item{A}{
The estimated rotation matrix.
}
\item{fitted}{
The fitted values in Euclidean coordinates (unit vectors).
}
}
\references{
Ted Chang (1986). Spherical Regression. Annals of Statistics, 14(3): 907--924.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{spher.cor}, \link{spml.reg}
}
}
\examples{
mx <- rnorm(3)
mx <- mx/sqrt( sum(mx^2) )
my <- rnorm(3)
my <- my/sqrt( sum(my^2) )
x <- rvmf(100, mx, 15)
A <- rotation(mx, my)
y <- x \%*\% t(A)
mod <- spher.reg(y, x, euclidean = TRUE)
A
mod$A ## exact match, no noise
y <- x \%*\% t(A)
y <- y + rvmf(100, colMeans(y), 40)
mod <- spher.reg(y, x, euclidean = TRUE)
A
mod$A ## noise added, more relistic example
}
\keyword{ Spherical-Spherical regression }
\keyword{ Spherical data }
