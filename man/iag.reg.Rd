\name{Spherical regression using rotationally symmetric distributions}
\alias{iag.reg}
\alias{vmf.reg}
\alias{sipc.reg}
\title{
Spherical regression using rotationally symmetric distributions
}

\description{
Spherical regression using rotationally symmetric distributions.
}

\usage{
iag.reg(y, x, con = TRUE, xnew = NULL, tol = 1e-06)
vmf.reg(y, x, con = TRUE, xnew = NULL, tol = 1e-06)
sipc.reg(y, x, con = TRUE, xnew = NULL, tol = 1e-06)
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
\item{tol}{
A tolerance value to decide when to stop the successive optimaizations.
}
}

\details{
The second parametrization of the projected normal and of the von Mises-Fisher regression (Paine et al., 2020) is applied. The same is true for the SIPC distribution. For more information see the paper by Paine et al. (2020).
}

\value{
A list including:
\item{loglik}{
The log-likelihood of the regression model.
}
\item{fit}{
This is a measure of fit of the estimated values, defined as \eqn{\sum_{i=1}^ny_i^T\hat{y}_i}. This appears if the argument "xnew" is NULL.
}
\item{beta}{
The beta coefficients.
}
\item{seb}{
The standard error of the beta coefficients.
}
\item{ki}{
The norm of the fitted values. In the von Mises-Fisher regression this is the concentration parameter of each observation. In the projected normal this are the norms of the fitted values before being projected onto the sphere. This is returned if the argument "xnew" is NULL.
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


Tsagris M. and Alzeley O. (2023). Circular and spherical projected Cauchy distributions:
A Novel Framework for Circular and Directional Data Modeling.
https://arxiv.org/pdf/2302.02468.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{esag.reg}, \link{vmfreg}, \link{spml.reg}
}
}

\examples{
y <- rvmf(150, rnorm(3), 5)
a1 <- iag.reg(y, iris[, 4])
a2 <- iag.reg(y, iris[, 4:5])

b1 <- vmf.reg(y, iris[, 4])
b2 <- vmf.reg(y, iris[, 4:5])
}

