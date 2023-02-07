\name{Logarithm of the Kent distribution normalizing constant}
\alias{kent.logcon}
\title{
Logarithm of the Kent distribution normalizing constant
}
\description{
Logarithm of the Kent distribution normalizing constant.
}

\usage{
kent.logcon(k, b, j = 100)
}

\arguments{
 \item{k}{
The conencration parameter, \eqn{\kappa}.
}
\item{b}{
The ovalness parameter, \eqn{\beta}.
}
\item{j}{
The number of the terms in the sum to use. By default this is 100.
}
}

\details{
It calculates logarithm of the normalising constant of the Kent distribution.
}

\value{
The value of the logarithm of the normalising constant of the Kent distribution.
}

\references{
Kent John (1982). The Fisher-Bingham distribution on the sphere. 
Journal of the Royal Statistical Society, Series B, 44(1): 71--80.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{fb.saddle}, \link{kent.mle}
}
}

\examples{
kent.logcon(10, 2)
fb.saddle( c(0, 10, 0), c(0, -2, 2) )
}

\keyword{ Kent distribution }
\keyword{ Normalising constant }
\keyword{ Saddlepoint approximation }

