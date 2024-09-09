\name{Spherical regression using the ESAG distribution}
\alias{esag.reg}
\title{
Spherical regression using the ESAG distribution
}

\description{
Spherical regression using the ESAG distribution.
}

\usage{
esag.reg(y, x, con = TRUE, xnew = NULL, lati = 10, longi = 10, tol = 1e-06)
}

\arguments{
\item{y}{
A matrix with 3 columns containing the (unit vector) spherical data.
}
\item{x}{
The predictor variable(s), they can be continnuous, spherical, categorical or a mix of them.
}
\item{con}{
Do you want the constant term in the regression?
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{lati}{
A positive number determing the range of degrees to move left and right from the latitude center. This number and the next determine the grid of points to search for the Q matrix described in Paine et al. (2020).
}
\item{longi}{
A positive number determing the range of degrees to move up and down from the longitude center. This number and the previous determine the grid of points to search for the Q matrix described in Paine et al. (2020).
}
\item{tol}{
A tolerance value to decide when to stop the successive optimizations.
}
}

\details{
The second parametrization of the ESAG regression (Paine et al., 2020) is applied.
}

\value{
A list including:
\item{loglik}{
The log-likelihood of the regression model.
}
\item{param}{
A vector with three numbers. A measure of fit of the estimated values, defined as \eqn{\sum_{i=1}^ny_i^T\hat{y}_i}. This appears if the argument "xnew" is NULL. The \eqn{\rho \in (0,1]} (smallest eigenvalue of the covariance matrix)), and the angle of rotation \eqn{psi}.
}
\item{gam}{
The two \eqn{\gamma} parameters.
}
\item{beta}{
The beta coefficients.
}
\item{seb}{
The standard error of the beta coefficients.
}
\item{est}{
The fitted values of xnew if "xnew" is NULL. If it is not NULL, the fitted values for the "xnew" you supplied will be returned.
}
}

\references{
P. J. Paine, S. P. Preston, M. Tsagris and Andrew T. A. Wood (2020).
Spherical regression models with general covariates and anisotropic errors.
Statistics and Computing, 30(1): 153--165.
https://link.springer.com/content/pdf/10.1007%2Fs11222-019-09872-2.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{esag.mle}, \link{iag.reg}, \link{spml.reg}
}
}

\examples{
y <- resag( 25, rnorm(3), c(1, 1) )
## this is a small example to pass CRAN's check because the default argument values
## of lati and longi require many seconds
a <- esag.reg(y, iris[1:25, 4], lati = 2, longi = 2)
}

