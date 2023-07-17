\name{Spherical and hyper-spherical distance correlation}
\alias{spher.dcor}

\title{Spherical and hyper-spherical distance correlation
}

\description{
Spherical and hyper-spherical distance correlation.
}

\usage{
spher.dcor(x, y)
}

\arguments{
\item{x}{
A matrix with directional data, i.e. unit vectors.
}
\item{y}{
A matrix with directional data, i.e. unit vectors.
}
}

\details{
The distance correlation between two spherical or hyper-spherical variables is computed.
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
\code{\link{circ.dcor}
}
}

\examples{
y <- rvmf(50, rnorm(3), 4)
x <- rvmf(50, rnorm(3), 4)
spher.dcor(x, y)
}

