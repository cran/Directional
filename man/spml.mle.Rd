\name{MLE of some circular distributions}
\alias{spml.mle}
\alias{wrapcauchy.mle}
\alias{circexp.mle}
\alias{circbeta.mle}
\alias{cardio.mle}
\alias{ggvm.mle}
\alias{cipc.mle}
\alias{gcpc.mle}
\title{
MLE of some circular distributions
}

\description{
MLE of some circular distributions.
}

\usage{
spml.mle(x, rads = FALSE, tol = 1e-07)
wrapcauchy.mle(x, rads = FALSE, tol = 1e-07)
circexp.mle(x, rads = FALSE, tol = 1e-06)
circbeta.mle(x, rads = FALSE)
cardio.mle(x, rads = FALSE)
ggvm.mle(phi, rads = FALSE)
cipc.mle(x, rads = FALSE)
gcpc.mle(x, rads = FALSE)
}

\arguments{
\item{x}{
A numerical vector with the circular data. They can either be expressed in radians or in degrees.
}
\item{phi}{
A numerical vector with the circular data. They can either be expressed in radians or in degrees.
}
\item{rads}{
If the data are in radians set this to TRUE.
}
\item{tol}{
The tolerance level to stop the iterative process of finding the MLEs.
}
}

\details{
The parameters of the bivariate angular Gaussian (spml.mle), wrapped Cauchy, circular exponential, cardioid, circular beta, 
geometrically generalised von Mises, CIPC (reparametrised version of the wrapped Cauchy) and GCPC (generalisation of the CIPC) distributions are estimated. For the Wrapped Cauchy, the iterative procedure described by Kent and Tyler (1988) is used. The Newton-Raphson algortihm for the angular Gaussian is described in the regression setting in Presnell et al. (1998). 
The circular exponential is also known as wrapped exponential distribution. 
}

\value{
A list including:
\item{iters}{
The iterations required until convergence.
}
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{param}{
A vector consisting of the estimates of the two parameters, the mean direction for both distributions
and the concentration parameter \eqn{kappa} and the \eqn{rho} for the von Mises and wrapped Cauchy respectively.
For the circular beta this contains the mean angle and the \eqn{\alpha} and \eqn{\beta} parameters. For the cardioid
distribution this contains the \eqn{\mu} and \eqn{rho} parameters. For the generalised von Mises this is a vector
consisting of the \eqn{\zeta}, \eqn{\kappa}, \eqn{\mu} and \eqn{\alpha} parameters of the generalised von Mises
distribution as described in Equation (2.7) of Dietrich and Richter (2017).
}
\item{gamma}{
The norm of the mean vector of the angular Gaussian, the CIPC and the GCPC distributions.
}
\item{mu}{
The mean vector of the angular Gaussian, the CIPC and the GCPC distributions.
}
\item{mumu}{
In the case of "angular Gaussian distribution this is the mean angle in radians.
}
\item{circmu}{
In the case of the CIPC and the GCPC this is the mean angle in radians.
}
\item{rho}{
For the GCPC distribution this is the eigenvalue of the covariance matrix, or the covariance determinant.
}
\item{lambda}{
The lambda parameter of the circular exponential distribution.
}
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra S. (2012). A short note on parameter approximation for von Mises-Fisher distributions:
and a fast implementation of \eqn{I_s(x)}. Computational Statistics, 27(1): 177--190.

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data.
Journal of the American Statistical Association, 93(443): 1068--1077.

Kent J. and Tyler D. (1988). Maximum likelihood estimation for the wrapped Cauchy distribution.
Journal of Applied Statistics, 15(2): 247--254.

Dietrich T. and Richter W. D. (2017). Classes of geometrically generalized von Mises distributions.
Sankhya B, 79(1): 21--59.

https://en.wikipedia.org/wiki/Wrapped_exponential_distribution

Jammalamadaka S. R. and Kozubowski T. J. (2003). A new family of circular models: 
The wrapped Laplace distributions.
Advances and Applications in Statistics, 3(1), 77-103.

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
\code{\link{circ.summary}, \link{purka.mle}, \link{rvonmises}, \link{vmf.mle}, \link{rvmf}
}
}

\examples{
x <- rvonmises(1000, 3, 9)
spml.mle(x, rads = TRUE)
wrapcauchy.mle(x, rads = TRUE)
circexp.mle(x, rads = TRUE)
ggvm.mle(x, rads = TRUE)
}

\keyword{ bivariate angular Gaussian }
\keyword{ wrapped Cauchy distribution }
\keyword{ circular data }
