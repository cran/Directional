\name{Projections based test of uniformity}
\alias{ptest}
\title{
Projections based test of uniformity
}

\description{
It checkes whether the data are uniformly distributed on the circle or the (hyper-)sphere.
}

\usage{
ptest(x, B = 100)
}

\arguments{
\item{x}{
A matrix containing the data, unit vectors.
}
\item{B}{
The number of random uniform projections to use.
}
}

\details{
For more details see Cuesta-Albertos, Cuevas and Fraiman (2009).
}

\value{
A list including:
\item{pvalues}{
The p-values of the Kolmogorov-Smirnov tests.
}
\item{pvalue}{
The p-value of the test based on the Benjamini and Heller (2008) procedure.
}
}

\references{
Cuesta-Albertos J. A., Cuevas A. and Fraiman, R. (2009).
On projection-based tests for directional and compositional data.
Statistics and Computing, 19: 367--380.

Benjamini Y. and Heller R. (2008). Screening for partial conjunction hypotheses.
Biometrics, 64(4): 1215--1222.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rayleigh}, \link{kuiper}
}
}

\examples{
x <- rvmf(100, rnorm(5), 1)  ## Fisher distribution with low concentration
ptest(x)
}
