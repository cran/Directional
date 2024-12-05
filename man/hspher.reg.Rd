\name{Hyper spherical-spherical regression}
\alias{hspher.reg}
\title{
Hyper spherical-spherical regression
}

\description{
Regression when both the dependent and independent variables are directional data-.
}

\usage{
hspher.reg(y, x, xnew = NULL)
}

\arguments{
\item{y}{
The dependent variable; a matrix with either two columns, latitude and longitude, either in radians or in degrees. Alternatively it is a matrix with three columns, unit vectors.
}
\item{x}{
The dependent variable; a matrix with either two columns, latitude and longitude, either in radians or in degrees. Alternatively it is a matrix with three columns, unit vectors. The two matrices must agree in the scale and dimensions.
}
\item{xnew}{
The new values of some directional independent variable(s) whose directional response values you want to predict. If you have no new x values, leave it NULL (default).
}
}

\details{
Spherical regression as proposed by Chang (1986) is implemented. If the estimated rotation matrix has a determinant equal to -1, singular value decomposition is performed and the last unit vector is multiplied by -1.
}

\value{
A list including:
\item{A}{
The estimated rotation matrix.
}
\item{est}{
The fitted values in unit vectors, if the argument xnew is not NULL.
}
}

\references{
Ted Chang (1986). Spherical Regression. Annals of Statistics, 14(3): 907--924.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{spher.cor}, \link{spml.reg}, \link{spher.reg}, \link{sphereplot}
}
}

\examples{
mx <- rnorm(5)
mx <- mx/sqrt( sum(mx^2) )
my <- rnorm(5)
my <- my/sqrt( sum(my^2) )
x <- rvmf(100, mx, 15)
A <- rotation(mx, my)
y <- x \%*\% t(A)
mod <- hspher.reg(y, x)
A
mod$A ## exact match, no noise
y <- x \%*\% t(A)
y <- y + rvmf(100, colMeans(y), 40)
mod <- hspher.reg(y, x)
A
mod$A ## noise added, more relistic example
}

