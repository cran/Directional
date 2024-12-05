\name{Score test for many simple CIPC and SMPL regressions}
\alias{score.cipc}
\alias{score.spml}
\title{
Score test for many simple CIPC and SPML regressions
}

\description{
Score test for many simple CIPC and SPML regressions.
}

\usage{
score.cipc(y, X, rads = TRUE, tol = 1e-06)
score.spml(y, X, rads = TRUE, tol = 1e-06)
}

\arguments{
\item{y}{
The dependent variable,  a numerical vector, it can be in radians or degrees.
}
\item{X}{
A matrix with many numerical independent variables.
}
\item{rads}{
If the dependent variable is expressed in rads, this should be TRUE and FALSE otherwise.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm in the null model (CIPC without covariates).
}
}

\details{
The score test uses the first derivative (score function) of the regression log-likelihood and it is asymptotically correct. So, this function requires sample sizes or at least 1,000 observations. The CIPC is basically the Wrapped Cauchy distribution (Tsagris and Alzeley, 2024) and SPML is the bivariate projected normal (Presnell et al., 1998).
}

\value{
A matrix with two columns, the test statistic and its associated p-value.
}

\references{
Tsagris M. and Alzeley O. (2024). Circular and spherical projected Cauchy distributions: A Novel Framework for Circular and Directional Data Modeling. Australian & New Zealand Journal of Statistics (accepted for publication). https://arxiv.org/pdf/2302.02468.pdf

Presnell B., Morrison S. P. and Littell R. C. (1998). Projected multivariate linear models for
directional data. Journal of the American Statistical Association, 93(443): 1068--1077.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cipc.reg}, \link{spml.reg}, \link{cipc.mle}, \link{spml.mle},
}
}

\examples{
y <- rcipc(500, omega = 2, g = 5)
x <- matrix( rnorm(500 * 10), ncol = 10 )
a <- score.cipc(y, x)
}


