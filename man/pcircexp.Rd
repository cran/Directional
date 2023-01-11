\name{Probability density function of the circular (wrapped) exponential distribution}
\alias{pcircexp}
\title{
Probability density function of the circular (wrapped) exponential distribution
}

\description{
Probability density function of the circular (wrapped) exponential distribution.
}

\usage{
pcircexp(u, lambda, rads = FALSE)
}

\arguments{
\item{u}{
A numerical value, either in radians or in degrees.
}
\item{lambda}{
The \eqn{\lambda} parameter. This must be positive.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
}

\details{
This value calculates the probability of x being less than u.
}

\value{
The probability that of x being less than u, where u follows the wrapped exponential distribution.
}

\references{
Jammalamadaka, S. R. & Kozubowski, T. J. (2003). A new family of circular models: The wrapped Laplace distributions.
Advances and Applications in Statistics, 3(1), 77-103.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{pvm}, \link{circ.summary}, \link{dvm}
}
}

\examples{
pcircexp(c(1, 2), 2, rads = TRUE)
}
