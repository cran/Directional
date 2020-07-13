\name{MLE of some circular distributions}
\alias{spml.mle}
\alias{wrapcauchy}
\alias{circexp.mle}
\title{
MLE of some circular distributions
}

\description{
MLE of some circular distributions.
}

\usage{
spml.mle(x, rads = FALSE, tol = 1e-07)
wrapcauchy(x, rads = FALSE, tol = 1e-07)
circexp.mle(x, rads = FALSE, tol = 1e-06)
}

\arguments{
\item{x}{
A numerical vector with the circular data. They must be expressed in radians.
}
\item{rads}{
If the data are in radians set this to TRUE.
}
\item{tol}{
The tolerance level to stop the iterative process of finding the MLEs.
}
}

\details{
The parameters of the bivariate angular Gaussian and wrapped Cauchy distributions are estimated. For the Wrapped Cauchy,
the iterative procedure described by Kent and Tyler (1988) is used. The Newton-Raphson algortihm for the angular Gaussian is described
in the regression setting in Presnell et al. (1998).
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
and the concentration parameter kappa and the rho for the von Mises and wrapped Cauchy respectively.
}
\item{gamma}{
The norm of the mean vector of the angular Gaussian distribution.
}
\item{mu}{
The mean vector of the angular Gaussian distribution.
}
\item{mumu}{
In the case of "spml.mle" this is the mean angle in radians.
}
\item{lambda}{
The lambda parameter of the circular exponential distribution.
}
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley \& Sons.

Sra S. (2012). A short note on parameter approximation for von Mises-Fisher distributions:
and a fast implementation of Is(x). Computational Statistics, 27(1): 177-190.

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data.
Journal of the American Statistical Association, 93(443): 1068-1077.

Kent J. and Tyler D. (1988). Maximum likelihood estimation for the wrapped Cauchy distribution.
Journal of Applied Statistics, 15(2): 247--254.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{vmf}, \link{rvonmises}, \link{rvmf}
}
}

\examples{
x <- rvonmises(1000, 3, 9)
spml.mle(x, rads = TRUE)
wrapcauchy(x, rads = TRUE)
circexp.mle(x, rads = TRUE)
}

\keyword{ bivariate angular Gaussian }
\keyword{ wrapped Cauchy distribution }
\keyword{ circular data }
