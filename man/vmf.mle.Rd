\name{MLE of (hyper-)spherical distributions}
\alias{vmf.mle}
\alias{multivmf.mle}
\alias{acg.mle}
\alias{iag.mle}
\alias{spcauchy.mle}
\title{
MLE of (hyper-)spherical distributions
}

\description{
MLE of (hyper-)spherical distributions.
}

\usage{
vmf.mle(x, fast = FALSE, tol = 1e-07)
multivmf.mle(x, ina, tol = 1e-07, ell = FALSE)
acg.mle(x, tol = 1e-07)
iag.mle(x, tol = 1e-07)
spcauchy.mle(x, tol = 1e-06)
}

\arguments{
\item{x}{
A matrix with directional data, i.e. unit vectors.
}
\item{fast}{
IF you want a faster version, but with fewer information returned, set this equal to TRUE.
}
\item{ina}{
A numerical vector with discrete numbers starting from 1, i.e. 1, 2, 3, 4,... or a factor variable. Each number denotes a sample or group.
If you supply a continuous valued vector the function will obviously provide wrong results.
}
\item{ell}{
This is for the multivmf.mle only. Do you want the log-likelihood returned? The default value is TRUE.
}
\item{tol}{
The tolerance value at which to terminate the iterations.
}
}

\details{
The vmf.mle estimates the mean direction and concentration of a fitted von Mises-Fisher distribution.

The von Mises-Fisher distribution for groups of data is also implemented.

The acg.mle fits the angular central Gaussian distribution. There is a constraint on the estimated covariance matrix; its trace is equal
to the number of variables. An iterative algorithm takes place and convergence is guaranteed.

The iag.mle implements MLE of the spherical projected normal distribution, for spherical data only.

The spcauchy.mle estimates the parameters of the spherical Cauchy distribution, for any dimension. The name sounds confusing, but it is implemented for
arbitrary dimensions, not only the sphere.
}

\value{
For the von Mises-Fisher a list including:
\item{loglik}{
The maximum log-likelihood value.
}
\item{mu}{
The mean direction.
}
\item{kappa}{
The concentration parameter.
}
For the multi von Mises-Fisher a list including:
\item{loglik}{
A vector with the maximum log-likelihood values if ell is set to TRUE. Otherwise NULL is returned.
}
\item{mi}{
A matrix with the group mean directions.
}
\item{ki}{
A vector with the group concentration parameters.
}
For the angular central Gaussian a list including:
\item{iter}{
The number if iterations required by the algorithm to converge to the solution.
}
\item{cova}{
The estimated covariance matrix.
}
For the spherical projected normal a list including:
\item{iters}{
The number of iteration required by the Newton-Raphson.
}
\item{mesi}{
A matrix with two rows. The first row is the mean direction and the second is the mean vector.
The first comes from the second by normalising to have unit length.
}
\item{param}{
A vector with the elements, the norm of mean vector, the log-likelihood and the log-likelihood of the spherical uniform distribution.
The third value helps in case you want to do a log-likelihood ratio test for uniformity.
}
For the spherical Cauchy a list including:
\item{mu}{
The mean direction.
}
\item{rho}{
The concetration parameter, this takes values in [0, 1).
}
\item{loglik}{
The log-likelihood value.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra, S. (2012). A short note on parameter approximation for von Mises-Fisher distributions: and a fast implementation of Is(x).
Computational Statistics, 27(1): 177--190.

Tyler D. E. (1987). Statistical analysis for the angular central Gaussian distribution on the sphere.
Biometrika 74(3): 579-589.

Paine P.J., Preston S.P., Tsagris M and Wood A.T.A. (2018). An Elliptically Symmetric Angular Gaussian Distribution.
Statistics and Computing, 28: 689-697.

Kato S. and McCullagh P. (2018). Mobius transformation and a Cauchy family on the sphere. arXiv preprint arXiv:1510.07679.
}

\author{
Michail Tsagris
R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{racg}, \link{vm.mle}, \link{rvmf}
}
}

\examples{
m <- c(0, 0, 0, 0)
s <- cov(iris[, 1:4])
x <- racg(100, s)
mod <- acg.mle(x)
mod
cov2cor(mod$cova)  ## estimated covariance matrix turned into a correlation matrix
cov2cor(s)  ## true covariance matrix turned into a correlation matrix
vmf.mle(x)
x <- rbind( rvmf(100,rnorm(4), 10), rvmf(100,rnorm(4), 20) )
a <- multivmf.mle(x, rep(1:2, each = 100) )
}

