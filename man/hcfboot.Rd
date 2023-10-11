\name{Bootstrap ANOVA for (hyper-)spherical data}
\alias{hcfboot}
\alias{hetboot}

\title{
Bootstrap ANOVA for (hyper-)spherical data
}

\description{
Bootstrap ANOVA for (hyper-)spherical data.
}

\usage{
hcfboot(x, ina, B = 999)
hetboot(x, ina, B = 999)
}

\arguments{
\item{x}{
A matrix with the combined data (from all groups) in Euclidean coordinates, i.e. unit vectors.
}
\item{ina}{
The grouping variables. A factor or a numerical vector specifying the groups to which each observation belongs to.
}
\item{B}{
The number of bootstraps to perform.
}
}

\details{
The high concentration (hcfboot), or the non equal
concentration parameters approach (hetboot) is used.
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

Tsagris M. and Alenazi A. (2022). An investigation of hypothesis testing procedures for circular and spherical mean vectors. Communications in Statistics-Simulation and Computation (Accepted for publication).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hcf.boot}, \link{hcf.aov}
}
}

\examples{
x <- rvmf(60, rnorm(3), 10)
ina <- rep(1:3, each = 20)
hcfboot(x, ina)
}
