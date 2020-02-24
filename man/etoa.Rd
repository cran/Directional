\name{Transform unit vectors to angular data}
\alias{etoa}

\title{
Transform unit vectors to angular data
}

\description{
Transform unit vectors to angular data.
}

\usage{
etoa(x)
}

\arguments{
\item{x}{
A numerical matrix with directional data, i.e. unit verctors.
}
}

\details{
from the Euclidean coordinates the data are mapped to angles, expressed in rads. 
}

\value{
A list including:
\item{mu}{
A matrix with angles. The number of columns is one less than that of the original matrix.
}
}

\references{
https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{vmnb.pred}, \link{weibull.nb} }
}

\examples{
x <- rvmf(10, rnorm(3), 5)
y <- etoa(x)
}
