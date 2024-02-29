\name{Density of some circular distributions}
\alias{dvm}
\alias{dspml}
\alias{dwrapcauchy}
\alias{dcircpurka}
\alias{dggvm}
\alias{dcircbeta}
\alias{dcardio}
\alias{dcircexp}
\alias{dcipc}
\alias{dgcpc}
\alias{dmmvm}
\title{
Density of some circular distributions
}

\description{
Density of some circular distributions.
}

\usage{
dvm(x, m, k, rads = FALSE, logden = FALSE)
dspml(x, mu, rads = FALSE, logden = FALSE)
dwrapcauchy(x, m, rho, rads = FALSE, logden = FALSE)
dcircpurka(x, m, a, rads = FALSE, logden = FALSE)
dggvm(x, param, rads = FALSE, logden = FALSE)
dcircbeta(x, m, a, b, rads = FALSE, logden = FALSE)
dcardio(x, m, rho, rads = FALSE, logden = FALSE)
dcircexp(x, lambda, rads = FALSE, logden = FALSE)
dcipc(x, omega, g, rads = FALSE, logden = FALSE)
dgcpc(x, omega, g, rho, rads = FALSE, logden = FALSE)
dmmvm(x, m, k, N, rads = FALSE, logden = FALSE)
}

\arguments{
\item{x}{
A vector with circular data.
}
\item{m}{
The mean value of the von Mises distribution and of the cardioid, a scalar.
This is the median for the circular Purkayastha distribution.
}
\item{mu}{
The mean vector, a vector with two values for the "spml".
}
\item{omega}{
The location parameter of the CIPC and GCPC distributions.
}
\item{g}{
The norm of the mean vector for the CIPC and GCPC distributions.
}
\item{k}{
The concentration parameter.
}
\item{rho}{
For the wrapped Cauchy and Cardioid distributions, this is the \eqn{\rho} parameter.
For the GCPC distribution this is the eigenvalue parameter, or covariance determinant.
}
\item{a}{
The \eqn{\alpha} parameter of the circular Purkayastha distribution or the \eqn{\alpha} parameter of the
circular Beta distribution.
}
\item{b}{
The \eqn{\beta} parameter of the circular Beta distribution.
}
\item{lambda}{
The \eqn{\lambda} parameter of the circular (or wrapped) exponential distribution. This must be positive.
}
\item{param}{
The vector of parameters of the GGVM distribution as returned by the function \code{\link{ggvm.mle}}.
}
\item{N}{
The number of modes to consider in the multi-modal von Mises distribution.
}
\item{rads}{
If the data are in rads, then this should be TRUE, otherwise FALSE.
}
\item{logden}{
If you the logarithm of the density values set this to TRUE.
}
}

\details{
The density of the von Mises, bivariate projected normal, cardio, circular exponential,
wrapped Cauchy, circular Purkayastha, CIPC or GCPC distributions is computed.
}

\value{
A vector with the (log) density values of x.
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Tsagris M. and Alzeley O. (2023). Circular and spherical projected Cauchy distributions:
A Novel Framework for Circular and Directional Data Modeling.
https://arxiv.org/pdf/2302.02468.pdf

Presnell B., Morrison S. P. and Littell R. C. (1998). Projected multivariate linear models for
directional data. Journal of the American Statistical Association, 93(443): 1068--1077.

Jammalamadaka S. R. and Kozubowski T. J. (2003). A new family of circular models:
The wrapped Laplace distributions. Advances and Applications in Statistics, 3(1): 77--103.

Barnett M. J. and Kingston R. L. (2024). A note on the Hendrickson-Lattman phase probability
distribution and its equivalence to the generalized von Mises distribution.
Journal of Applied Crystallography, 57(2).

Paula F. V., Nascimento A. D., Amaral G. J. and Cordeiro G. M. (2021).
Generalized Cardioid distributions for circular data analysis. Stats, 4(3): 634--649.

Zheng Sun (2009). Comparing measures of fit for circular distributions. MSc Thesis, University of Victoria.
file:///C:/Users/mtsag/Downloads/zhengsun_master_thesis.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{dkent}, \link{rvonmises}, \link{desag}
}
}

\examples{
x <- rvonmises(500, m = 2.5, k = 10, rads = TRUE)
mod <- circ.summary(x, rads = TRUE, plot = FALSE)
den <- dvm(x, mod$mesos, mod$kappa, rads = TRUE, logden = TRUE )
mod$loglik
sum(den)
}

