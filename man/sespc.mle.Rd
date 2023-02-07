\name{MLE of the SESPC distribution}
\alias{sespc.mle}
\title{
MLE of the SESPC distribution
}

\description{
MLE of the SESPC distribution.
}

\usage{
sespc.mle(y, full = FALSE, tol = 1e-06)
}

\arguments{
\item{y}{
A matrix with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{full}{
If you want some extra information, the inverse of the covariance matrix, set this equal to TRUE. Otherwise leave it FALSE.
}
\item{tol}{
A tolerance value to stop performing successive optimizations.
}
}

\details{
MLE of the SESPC distribution is implemented. SESPC stands for Spherical Elliptically Symmetric Projected Cauchy
and it was suugested by Tsagris and Alzeley (2023). Unlike the spherical independent projected Cauchy distribution this is rotationally symmetric and is a competitor of the spherical ESAG and Kent distributions (which are also ellitpically symmetric).
}

\value{
A list including:
\item{mu}{
The mean vector in \eqn{R^3}.
}
\item{theta}{
The two \eqn{\theta} parameters.
}
\item{loglik}{
The log-likelihood value.
}
\item{vinv}{
The inverse of the covariance matrix. It is returned if the argument "full" is TRUE.
}
\item{lambda}{
The \eqn{\lambda_2} parameter (smallest eigenvalue of the covariance matrix). It is returned if the argument "full" is TRUE.
}
\item{psi}{
The angle of rotation \eqn{\psi} set this equal to TRUE. It is returned if the argument "full" is TRUE.
}
\item{sipc.loglik}{
The log-likelihood value of the isotropic prohected Cuchy distribution, which is rotational symmetric.
}
}

\references{
Tsagris M. and Alzeley O. (2023). Circular and spherical projected Cauchy distributions.
https://arxiv.org/pdf/2302.02468.pdf

Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dsespc}, \link{rsespc}, \link{sipc.mle}, \link{esag.mle}, \link{spher.sespc.contour}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- rsespc(1000, m, c(1,0.5) )
sespc.mle(y)
}

\keyword{ SESPC distribution }
\keyword{ maximum likelihood estimation }
\keyword{ spherical data }
