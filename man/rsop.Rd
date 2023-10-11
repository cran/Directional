\name{Random sample of matrices in SO(p)}
\alias{rsop}
\title{
Random sample of matrices in SO(p)
}

\description{
Random sample of matrices in SO(p).
}

\usage{
rsop(n, p)
}

\arguments{
\item{n}{
The sample size, the number of matrices you want to generate.
}
\item{p}{
The dimensionality of the matrices.
}
}

\details{
The idea is very simple. Start with a unit vector pointing at the north pole (1,0,...,0). Then generate random numbers from a standard normal and scale them so that they have a unit length. To put it differently, a sample of n values from the uniform distribution on the sphere is generated. Then calculate the rotation matrix required to go from the north pole to each of a generated vector.
}

\value{
If n = 1 one matrix is returned. If n is greater than 1, an array with n matrices inside.
}

\references{
Amaral G.J.A., Dryden I.L. and Wood A.T.A. (2007).
Pivotal Bootstrap Methods for k-Sample Problems in Directional Statistics and Shape Analysis. 
Journal of the American Statistical Association, 102(478): 695--707.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rotation}, \link{Arotation}, \link{rot.matrix}
}
}

\examples{
x1 <- rsop(1, 3)
x2 <- rsop(10, 3)
x3 <- rsop(100, 10)
}

\keyword{ SO(p) }
\keyword{ Rotation matrix }
