\name{MLE of the ESAG distribution}
\alias{ESAGmle}
\title{
MLE of the ESAG distribution
}

\description{
MLE of the ESAG distribution.
}

\usage{
ESAGmle(y, full = FALSE, tol = 1e-06)
}

\arguments{
\item{y}{
A matrix with the data expressed in Euclidean coordinates, i.e. unit vectors.
}
\item{full}{
If you want some extra information, the inverse of the covariance matrix, the \eqn{rho} parameter (smallest eigenvalue of the covariance matrix) and the angle of rotation \eqn{\psi} set this equal to TRUE. Otherwise leave it FALSE.
}
\item{tol}{
A tolerance value to stop performing successive optimizations.
}
}

\details{
MLE of the MLE of the ESAG distributiontribution, on the sphere, is implemented. ESAG stands for Elliptically Symmetric Angular Gaussian
and it was suugested by Paine et al. (2017). Unlike the projected normal distribution this is rotationally symmetric and is a competitor
of the spherical Kent distribution (which is also non rotational symmetric).
}

\value{
A list including:
\item{mu}{
The mean vector in \eqn{R^3}.
}
\item{gam}{
The two gamma parameters.
}
\item{loglik}{
The log-likelihood value.
}
\item{vinv}{
The inverse of the covariance matrix. It is returned if the argument "full" is TRUE.
}
\item{rho}{
The \eqn{rho} parameter (smallest eigenvalue of the covariance matrix). It is returned if the argument "full" is TRUE.
}
\item{psi}{
The angle of rotation \eqn{\psi} set this equal to TRUE. It is returned if the argument "full" is TRUE.
}
\item{iag.loglik}{
The log-likelihood value of the isotropic angular Gaussian distribution. That is, the projected
normal distribution which is rotational symmetric.
}
}

\references{
Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.

Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ESAGdensity}, \link{ESAGsim}, \link{iag.mle}, \link{kent.mle}, \link{acg}, \link{circ.summary}, \link{sphereplot}
}
}

\examples{
m <- colMeans( as.matrix( iris[,1:3] ) )
y <- ESAGsim(1000, c(m, 1,0.5) )
ESAGmle(y)
}

\keyword{ ESAG distribution }
\keyword{ maximum likelihood estimation }
\keyword{ spherical data }
