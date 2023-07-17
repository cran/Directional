\name{Circular distance correlation between two circular variables}
\alias{circ.dcor}

\title{Circular distance correlation between two circular variables
}

\description{
Circular distance correlation between two circular variables.
}

\usage{
circ.dcor(theta, phi, rads = FALSE)
}

\arguments{
\item{theta}{
The first cirular variable.
}
\item{phi}{
The other cirular variable.
}
\item{rads}{
If the data are expressed in rads, then this should be TRUE. If the data are in degrees, then this is FALSE.
}
}

\details{
The angular data are transformed to their Euclidean coordinates and then the distance correlation is computed.
}

\value{
A list including:
\item{dcov}{
The distance covariance.
}
\item{dvarX}{
The distance variance of x.
}
\item{dvarY}{
The distance variance of Y.
}
\item{dcor}{
The distance correlation.
}
}

\references{
G.J. Szekely, M.L. Rizzo and N. K. Bakirov (2007). Measuring and Testing Independence 
by Correlation of Distances. Annals of Statistics, 35(6):2769-2794.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{circlin.cor}, \link{circ.cor2}, \link{spher.dcor}
}
}

\examples{
y <- runif(50, 0, 2 * pi)
x <- runif(50, 0, 2 * pi)
circ.dcor(x, y, rads = TRUE)
}

