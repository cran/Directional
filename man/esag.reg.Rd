\name{Spherical regression using the ESAG distribution}
\alias{esag.reg}
\title{
Spherical regression using the ESAG distribution
}

\description{
Spherical regression using the ESAG distribution.
}

\usage{
esag.reg(y, x, con = TRUE, xnew = NULL, n.axes = 50, Q = NULL)
}

\arguments{
\item{y}{
A matrix with 3 columns containing the (unit vector) spherical data.
}
\item{x}{
The predictor variable(s), they can be continuous, spherical, categorical or a mix of them.
}
\item{con}{
Do you want the constant term in the regression?
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
\item{n.axes}{
The approximate number of axes in the grid searched over for the matrix \eqn{\bm{Q}} described in Paine et al. (2020). The axes are evenly spread over a hemisphere. Larger values give a finer search, and the computational time grows linearly with this number.
}
\item{Q}{
If you already have a \eqn{3 \times 3} rotation matrix \eqn{\bm{Q}}, for example from a fit of a larger model, supply it here and no search is performed. This is useful when fitting nested models for likelihood ratio tests, since the same \eqn{\bm{Q}} should be used for all of them. If NULL (default), \eqn{\bm{Q}} is chosen by the grid search.
}
}

\details{
The second model structure of Paine et al. (2020) is applied:
\eqn{\bm{Q}^T\bm{y}_i \sim ESAG\left(\bm{\mu}(\bm{x}_i), \bm{\gamma}\right)}, with \eqn{\bm{\mu}(\bm{x}_i)=\bm{B}^T\bm{x}_i}, and the anisotropy parameters \eqn{\bm{\gamma}} constant, i.e. they do not depend on the covariates.

The \eqn{\gamma} parameters are defined relative to a frame built from the first coordinate axis, \eqn{\bm{e}_1}. Hence, without the rotation \eqn{\bm{Q}}, the anisotropic model would depend on the coordinate system of the data. \eqn{\bm{Q}} is treated as a tuning parameter and is chosen to maximise the log-likelihood. The log-likelihood depends on \eqn{\bm{Q}} only through the axis \eqn{\pm \bm{Q}\bm{e}_1}, so the search is performed over a grid of axes on a hemisphere (a Fibonacci grid), and for each axis the model is fitted by maximum likelihood.

All returned quantities, except for \eqn{\gamma}, \eqn{\rho} and \eqn{\psi}, are expressed in the original coordinate system of the data. The \eqn{\gamma} parameters, and hence \eqn{\rho} and \eqn{\psi}, are relative to the frame induced by the selected \eqn{\bm{Q}}.

If the errors are isotropic (\eqn{\gamma=0}), the matrix \eqn{\bm{Q}} plays no role. In that case use \code{\link{iag.reg}} instead.
}

\value{
A list including:
\item{loglik}{
The log-likelihood of the final regression model, i.e. the one fitted with the selected \eqn{\bm{Q}}.
}
\item{fit}{
A measure of fit of the estimated values, defined as \eqn{\sum_{i=1}^n\bm{y}_i^T\hat{\bm{y}}_i}.
}
\item{beta}{
The beta coefficients \eqn{\bm{B}} of the mean direction, in the original coordinate system.
}
\item{seb}{
The standard errors of the beta coefficients.
}
\item{gamma}{
The two \eqn{\gamma} parameters.
}
\item{seg}{
The standard errors of the \eqn{\gamma} parameters.
}
\item{rho}{
The \eqn{\rho \in (0,1]} value (the smallest eigenvalue of the covariance matrix).
}
\item{psi}{
The angle of rotation \eqn{\psi}.
}
\item{rl}{
The norms of the fitted \eqn{\bm{\mu}(\bm{x}_i)} vectors, i.e. the concentration of each observation.
}
\item{est}{
The predicted mean directions for "xnew". This is NULL if "xnew" is NULL.
}
\item{Q}{
The selected rotation matrix \eqn{\bm{Q}}.
}
\item{covb}{
The covariance matrix of all estimated parameters (first the \eqn{\gamma} parameters, then the beta coefficients).
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
\code{\link{esag.mle}, \link{sespc.reg}, \link{iag.reg}, \link{spml.reg}
}
}

\examples{
y <- resag( 25, rnorm(3), c(1, 1) )
## a small grid, so that the example runs quickly for CRAN's checks
a <- esag.reg(y, iris[1:25, 4], n.axes = 4)
}
