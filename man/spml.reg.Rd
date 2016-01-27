\name{spml.reg}
\alias{spml.reg}
\title{
Circular or angular regression.
}
\description{
Regression with circualr dependent variable and Euclidean or categorical independent variables.
}
\usage{
spml.reg(y, x, xnew = x, pred = F, rads = T)
}
\arguments{
  \item{y}{
  The dependent variable,  it can be in radians or degrees.
}
  \item{x}{
  The independent variable(s). Can be Euclidean or categorical. If they are categorical, the user must provide them as dummy variables. It does not accept factor variables.
}
  \item{xnew}{
  The new values of some independent variable(s) whose circular values you want to predict. Can be Euclidean or categorical. If they are categorical, the user must provide them as dummy variables. It does not accept factor variables.
}
  \item{pred}{
  A boolean variable taking values TRUE if the xnew is new data or FALSE if the xnew is the same as x. In the second case, the fitted values will be calculated.
}
  \item{rads}{
If the dependent variable is expressed in rads, this should be TRUE and FALSE otherwise.
}
}
\details{
The E-M algorithm as described in Presnell et al. (1998).
}
\value{
A list including:
\item{beta}{
The regression coefficients.
}
\item{seb}{
The standard errors of the coefficients.
}
\item{rho}{
The circular correlation of the ubeserved and the fitted data.
}
\item{R.sq}{
The squared circular correlation of the ubeserved and the fitted data.
}
\item{loglik}{
The vlaue of the maximised log-likelihood.
}
\item{fitted}{
The fitted values expressed in radians if the obsereved data are in radians and in degrees otherwise.
}

}
\references{
Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data. Journal of the American Statistical Association, 93(443): 1068-1077.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{circlin.cor}
}
}
\examples{
x <- rnorm(100)
z <- cbind(3+ 2 * x, 1 -3 * x)
y <- cbind( rnorm(100,z[ ,1], 1), rnorm(100, z[ ,2], 1) )
y <- y / sqrt( rowSums(y^2) )
y <- ( atan( y[, 2] / y[, 1] ) + pi * I(y[, 1] < 0) ) \%\% (2 * pi)
spml.reg(y, x, rads = TRUE)
}
\keyword{ Circular regression }
\keyword{ Projected normal }
