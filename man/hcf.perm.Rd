\name{Permutation based 2-sample mean test for (hyper-)spherical data}
\alias{hcf.perm}
\alias{lr.perm}
\alias{hclr.perm}
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
hclr.perm(x1, x2, B = 999)
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
The high concentration (hcf.perm), log-likelihood ratio (lr.perm), high concentration
log-likelihood ratio (hclr.perm), embedding approach (embed.perm) or the non equal
concentration parameters approach (het.perm) is used.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
The degrees of freedom of the test. Since these are permutation based tests this is "NA".
}
\item{p.value}{
The p-value of the test.
}
\item{alternative}{
A character with the alternative hypothesis.
}
\item{method}{
A character with the test used.
}
\item{data.name}{
A character vector with two elements.
}
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics.
Chicester: John Wiley & Sons.

Rumcheva P. and Presnell B. (2017). An improved test of equality of mean directions for the
Langevin-von Mises-Fisher distribution. Australian & New Zealand Journal of Statistics, 59(1), 119--135.

Tsagris M. and Alenazi A. (2022). An investigation of hypothesis testing procedures for circular
and spherical mean vectors. Communications in Statistics-Simulation and Computation (Accepted for publication).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hcf.boot}, \link{hcf.aov}, \link{spherconc.test}, \link{conc.test}
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
