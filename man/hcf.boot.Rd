\name{Bootstrap 2-sample mean test for (hyper-)spherical data}
\alias{hcf.boot}
\alias{lr.boot}
\alias{embed.boot}
\alias{het.boot}

\title{
Bootstrap 2-sample mean test for (hyper-)spherical data
}

\description{
Bootstrap 2-sample mean test for (hyper-)spherical data.
}

\usage{
hcf.boot(x1, x2, fc = TRUE, B = 999)

lr.boot(x1, x2, B = 999)

embed.boot(x1, x2, B = 999)

het.boot(x1, x2, B = 999)
}

\arguments{
\item{x1}{
A matrix with the data in Euclidean coordinates, i.e. unit vectors.
}
\item{x2}{
A matrix with the data in Euclidean coordinates, i.e. unit vectors.
}
\item{fc}{
A boolean that indicates whether a corrected F test should be used or not.
}
\item{B}{
The number of permutations to perform.
}
}

\details{
The high concentration (hcf.boot), log-likelihood ratio (lr.boot), embedding approach (embed.boot) or the non equal concentration parameters approach (het.boot) is used.
}

\value{
A vector including:
\item{test}{
The test statistic value.
}
\item{p-value}{
The p-value of the F test.
}
\item{kappa}{
The common concentration parameter kappa based on all the data.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hcf.aov}, \link{spherconc.test}, \link{conc.test}
}
}

\examples{
x <- rvmf(60, rnorm(3), 15)
ina <- rep(1:2, each = 30)
x1 <- x[ina == 1, ]
x2 <- x[ina == 2, ]
hcf.boot(x1, x2)
lr.boot(x1, x2)
het.boot(x1, x2)
}
