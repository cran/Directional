\name{Spherical regression using the SESPC distribution}
\alias{sespc.reg}
\title{
Spherical regression using the SESPC distribution
}

\description{
Spherical regression using the SESPC distribution.
}

\usage{
sespc.reg(y, x, con = TRUE, xnew = NULL, lati = 10, longi = 10, tol = 1e-06)
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
A positive number determing the range of degrees to move left and right from the latitude center. This number and the next determine the grid of points to search for the Q matrix described in Tsagris and Alzeley (2024).
}
\item{longi}{
A positive number determing the range of degrees to move up and down from the longitude center. This number and the previous determine the grid of points to search for the Q matrix described in Tsagris and Alzeley (2024).
}
\item{tol}{
A tolerance value to decide when to stop the successive optimizations.
}
}

\details{
Regression based on the SESPC distribution (Tsagris and Alzeley, 2024) is applied.
}

\value{
A list including:
\item{loglik}{
The log-likelihood of the regression model.
}
\item{param}{
A vector with three numbers. A measure of fit of the estimated values, defined as \eqn{\sum_{i=1}^ny_i^T\hat{y}_i}. This appears if the argument "xnew" is NULL. The \eqn{\rho \in (0,1]} (smallest eigenvalue of the covariance matrix)), and the angle of rotation \eqn{psi}.
}
\item{theta}{
The two \eqn{\theta} parameters.
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
Tsagris M. and Alzeley O. (2024). Circular and spherical projected Cauchy distributions: A Novel Framework for Circular and Directional Data Modeling. Australian & New Zealand Journal of Statistics (accepted for publication). https://arxiv.org/pdf/2302.02468.pdf
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
y <- rsespc( 150, rnorm(3), c(1, 1) )
## this is a small example to pass CRAN's check because the default argument values
## of lati and longi require many seconds
a <- sespc.reg(y, iris[, 4], lati = 2, longi = 2)
}

