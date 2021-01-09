\name{MLE of the generalised von Mises distribution}
\alias{ggvm.mle}
\title{
MLE of the geometrically generalised von Mises distribution
}

\description{
MLE of the geometrically generalised von Mises distribution.
}

\usage{
ggvm.mle(phi, rads = FALSE)
}

\arguments{
\item{phi}{
A numerical vector with the circular data.
}
\item{rads}{
Whether the data are in rads (TRUE) or not (FALSE).
}
}

\details{
The generalised von Mises distribution (Dietrich and Richter, 2016) is fitted to some data and its parameters are estimated.
}

\value{
A list including:
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{param}{
A vector consisting of the \eqn{\zeta}, \eqn{\kappa}, \eqn{\mu} and \eqn{\alpha} parameters of the generalised von Mises distribution as describe din Equation (2.7) of Dietrich and Richter (2016).
}
}

\references{
Dietrich, T. & Richter, W. D. (2017). Classes of geometrically generalized von Mises distributions. Sankhya B, 79(1): 21-59.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{circ.summary}, \link{rvonmises}
}
}

\examples{
x <- rvonmises(100, 2, 25, rads = TRUE)
circ.summary(x, rads = TRUE)
ggvm.mle(x, rads = TRUE)
}

\keyword{ Generalised von Mises distribution }
