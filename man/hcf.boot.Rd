\name{Bootstrap 2-sample mean test for (hyper-)spherical data}
\alias{hcf.boot}
\alias{lr.boot}
\alias{hclr.boot}
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
hclr.boot(x1, x2, B = 999)
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
The number of bootstraps to perform.
}
}

\details{
The high concentration (hcf.boot), log-likelihood ratio (lr.boot), high concentration
log-likelihood ratio (hclr.boot), embedding approach (embed.boot) or the non equal
concentration parameters approach (het.boot) is used.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
The degrees of freedom of the test. Since these are bootstrap based tests this is "NA".
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

Rumcheva P. and Presnell B. (2017). An improved test of equality of mean directions
for the Langevin-von Mises-Fisher distribution. Australian & New Zealand Journal of Statistics, 59(1): 119--135.

Tsagris M. and Alenazi A. (2024). An investigation of hypothesis testing procedures for circular
and spherical mean vectors. Communications in Statistics-Simulation and Computation, 53(3): 1387--1408.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hcf.aov}, \link{hcf.perm}, \link{hcfboot}
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
