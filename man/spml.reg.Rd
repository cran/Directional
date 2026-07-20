\name{Circular or angular regression}
\alias{spml.reg}
\alias{circpurka.reg}
\alias{cipc.reg}
\alias{gcpc.reg}
\title{
Circular or angular regression
}

\description{
Regression with circular dependent variable and Euclidean or categorical independent variables.
}

\usage{
spml.reg(y, x, rads = TRUE, xnew = NULL, seb = FALSE, tol = 1e-6, maxiters = 100)
circpurka.reg(y, x, rads = TRUE, xnew = NULL, tol = 1e-6, maxiters = 300)
cipc.reg(y, x, rads = TRUE, xnew = NULL, tol = 1e-6, maxiters = 100)
gcpc.reg(y, x, rads = TRUE, xnew = NULL)
}

\arguments{
\item{y}{
The dependent variable,  a numerical vector, it can be in radians or degrees, or a numerical
matrix with two columns containing eigenvectors.
}
\item{x}{
The independent variable(s). Can be Euclidean or categorical (factor variables).
}
\item{rads}{
If the dependent variable is expressed in rads, this should be TRUE and FALSE otherwise.
}
\item{xnew}{
The new values of some independent variable(s) whose circular values you want to predict.
Can be Euclidean or categorical. If you have no new x values, leave it NULL (default).
}
\item{seb}{
a boolean variable. If TRUE, the standard error of the coefficients will be be returned.
Set to FALSE in case of simulation studies or in other cases such as a forward regression
setting for example. In these cases, it can save some time.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations allowed in the Newton-Raphson algorithm.
}
}

\details{
For the spml.reg(), the Newton-Raphson algorithm is fitted in this regression as described in
Presnell et al. (1998). For the cipc.reg(), the Newton-Raphson algorithm is fitted in this regression
as described in Tsagris and Alenazy (2023). Note that the cipc.reg() is the same as the wrapped Cauchy
regression. For the circpurka.reg() a variant of the Newton-Raphson algorithm is used, followed by
the optim() function is employed. For the gcpc.reg() the optim() and the optimise() functions are
being used.
}

\value{
A list including:
\item{runtime}{
The runtime of the procedure.
}
\item{iters}{
The number of iterations required until convergence of the Newton-Raphson algorithm.
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
The fitted values expressed in radians if the obsereved data are in radians and in degrees otherwise.
If xnew is not NULL, i.e. if you have new x values, then the predicted values of y will be returned.
}
}

\references{
Presnell B., Morrison S. P. and Littell Ramon C. (1998). Projected multivariate linear models for
directional data.
Journal of the American Statistical Association, 93(443): 1068--1077.

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through
Maximum Likelihood Characterization.
The Indian Journal of Statistics, Series A, 53(1): 70--83

Tsagris M. and Alzeley O. (2025). Circular and spherical projected Cauchy distributions:
A Novel Framework for Circular and Directional Data Modeling.
Australian & New Zealand Journal of Statistics, 67(1): 77--103.
https://arxiv.org/pdf/2302.02468.pdf

Alzeley O. and Tsagris M. (2026). On the generalized circular projected Cauchy distribution.
Mathematics, 14(11): 1934.
https://www.mdpi.com/2227-7390/14/11/1934

}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{circlin.cor}, \link{circ.cor1}, \link{circ.cor2}, \link{spher.cor}, \link{spher.reg}
}
}

\examples{
x <- rnorm(100)
z <- cbind(3 + 2 * x, 1 -3 * x)
y <- cbind( rnorm(100,z[ ,1], 1), rnorm(100, z[ ,2], 1) )
y <- y / sqrt( rowSums(y^2) )
y <- ( atan( y[, 2] / y[, 1] ) + pi * I(y[, 1] < 0) ) \%\% (2 * pi)
a <- spml.reg(y, x, rads = TRUE, xnew = x)
b <- cipc.reg(y, x, rads = TRUE, xnew = x)
}


