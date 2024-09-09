\name{(Hyper-)spherical regression using rotational symmetric distributions}
\alias{vmfreg}
\alias{spcauchy.reg}
\alias{pkbd.reg}
\alias{pkbd.reg2}
\title{
(Hyper-)spherical regression using the rotational symmetric distributions
}

\description{
(Hyper-)spherical regression using the rotational symmetric distributions.
}

\usage{
vmfreg(y, x, con = TRUE, xnew = NULL, tol = 1e-06)
spcauchy.reg(y, x, con = TRUE, xnew = NULL, tol = 1e-06)
pkbd.reg(y, x, con = TRUE, xnew = NULL, tol = 1e-6)
pkbd.reg2(y, x, con = TRUE, xnew = NULL, tol = 1e-6)
}

\arguments{
\item{y}{
A matrix with any number of columns containing the (unit vector) (hyper-)spherical data.
}
\item{x}{
The predictor variable(s), they can be continnuous, (hyper-)spherical, categorical or a mix of them.
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
The second parametrization of the projected normal and of the von Mises-Fisher regression (Paine et al., 2020) is applied. The same is true for the SIPC distribution. For more information see the paper by Paine et al. (2020). The difference from \code{\link{vmf.reg}} is that the latter is designed for the sphere only, whereas this function works in the hyper-sphere also.

As for the spcauchy.reg() and pkbd.reg() they are based upon the spherical Cauchy (Kato and McCullagh, 2020) and the Poisson-kernel based (Golzy and Markatou, 2020) distributions. These two use Newton-Raphson, but the pkbd.reg2() uses the \code{\link{optim}}. We have noticed some numerical issues with the pkbd.reg() when the dimensionalities of the variables are large and this is why we also provide the (much slower) pkbd.reg2() function.
}

\value{
A list including:
\item{runtime}{
The runtime of the regression.
}
\item{iters}{
The number of iterations required until convergence of the Newton-Raphson algorithm.
}
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
The norm of the fitted values. In the von Mises-Fisher regression this is the concentration parameter of each observation. This is returned if the argument "xnew" is NULL.
}
\item{g2}{
The norm of the fitted values. In the spherical Cauchy and the PKBD regression this is the concentration parameter of each observation. This is returned if the argument "xnew" is NULL.
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

Kato S. and McCullagh P. (2020). Some properties of a Cauchy family on the sphere derived from the Mobius
transformations. Bernoulli, 26(4): 3224--3248.

Golzy M. and Markatou M. (2020). Poisson kernel-based clustering on the sphere: convergence properties,
identifiability, and a method of sampling. Journal of Computational and Graphical Statistics, 29(4): 758--770.

Tsagris M. (2024). Directional data analysis using the spherical Cauchy and the Poisson-kernel based distribution.
https://arxiv.org/pdf/2409.03292
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{ \link{esag.reg}, \link{vmf.reg}, \link{spml.reg}
}
}

\examples{
y <- rvmf(150, rnorm(5), 5)
a <- vmfreg(y, iris[, 1])
b <- spcauchy.reg(y, iris)
}

