\name{Unit vector(s) with a given angle}
\alias{vec}
\title{
Unit vector(s) with a given angle
}

\description{
Unit vector with a given angle from a given unit vector.
}

\usage{
vec(x, n = 1, deg = 90)
}

\arguments{
\item{x}{
A unit vector. If it is not a unit vector it becomes a unit vector.
}
\item{n}{
The number of unit vectors to return.
}
\item{deg}{
The angle between the given vector and the n vectors to be returned. This must be in degrees and it has to be between 0 and 180 degrees.
If the angle is 0, the same unit vector will be returned. If the angle is 180, the same unit vector with the signs changed will be returned.
}
}

\details{
The user provides a unit vector and the degrees. The function will return n unit vectors whose angle with the given unit vector equals the degrees given.
For example, if you want 10 unit vectors purpendicualr to the x put vec(x, 10, 90).
}

\value{
A list including:
\item{runtime}{
The runtime of the procedure.
}
\item{crit}{
The calculated angle between the given unit vector and each of the generated unit vectors.
}
\item{mat}{
A matrix with the n unit vectors.
}
}

\author{
Michail Tsagris
R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rvmf}, \link{rbingham}, \link{rfb}
}
}

\examples{
x <- rnorm(10)
a <- vec(x, 20, 90)
}

\keyword{ unit vectors }

