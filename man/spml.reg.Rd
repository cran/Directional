\name{spml.reg}
\alias{spml.reg}
\title{
Circular or angular regression.
}
\description{
Regression with circualr dependent variable and Euclidean or categorical independent variables.
}
\usage{
spml.reg(y, x, rads = TRUE, xnew = NULL)
}
\arguments{
  \item{y}{
  The dependent variable,  it can be in radians or degrees.
}
  \item{x}{
  The independent variable(s). Can be Euclidean or categorical (factor variables).
}
  \item{rads}{
If the dependent variable is expressed in rads, this should be TRUE and FALSE otherwise.
}
  \item{xnew}{
  The new values of some independent variable(s) whose circular values you want to predict. Can be Euclidean or categorical.
If they are categorical, the user must provide them as dummy variables. It does not accept factor variables. If you have no new x values, leave it NULL (default).
}
}
\details{
The E-M algorithm is fitted in this regression as described in Presnell et al. (1998).
}
\value{
A list including:
\item{runtime}{
The runtime of the procedure.
}
\item{beta}{
The regression coefficients.
}
\item{seb}{
The standard errors of the coefficients.
}
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{est}{
The fitted values expressed in radians if the obsereved data are in radians and in degrees otherwise. If xnew is not NULL, i.e. if you have new x values, then the predicted values of y will be returned.
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
