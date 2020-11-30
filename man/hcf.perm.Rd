\name{Permutation based 2-sample mean test for (hyper-)spherical data}
\alias{hcf.perm}
\alias{lr.perm}
\alias{embed.perm}
\alias{het.perm}

\title{
Permutation based 2-sample mean test for (hyper-)spherical data
}

\description{
Permutation based 2-sample mean test for (hyper-)spherical data.
}

\usage{
hcf.perm(x1, x2, B = 999)

lr.perm(x1, x2, B = 999)

embed.perm(x1, x2, B = 999)

het.perm(x1, x2, B = 999)
}

\arguments{
\item{x1}{
A matrix with the data in Euclidean coordinates, i.e. unit vectors.
}
\item{x2}{
A matrix with the data in Euclidean coordinates, i.e. unit vectors.
}
\item{B}{
The number of permutations to perform.
}
}

\details{
The high concentration (hcf.perm), log-likelihood ratio (lr.perm), embedding approach (embed.perm) or the non equal concentration parameters approach (het.perm) is used.
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
\code{\link{hcf.boot}, \link{spherconc.test}, \link{conc.test}
}
}

\examples{
x <- rvmf(60, rnorm(3), 15)
ina <- rep(1:2, each = 30)
x1 <- x[ina == 1, ]
x2 <- x[ina == 2, ]
hcf.perm(x1, x2)
lr.perm(x1, x2)
het.boot(x1, x2)
}
