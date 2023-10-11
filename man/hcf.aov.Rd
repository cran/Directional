\name{Anova for (hyper-)spherical data}
\alias{hcf.aov}
\alias{lr.aov}
\alias{embed.aov}
\alias{het.aov}
\alias{hclr.aov}

\title{
Analysis of variance for (hyper-)spherical data
}

\description{
Analysis of variance for (hyper-)spherical data.
}

\usage{
hcf.aov(x, ina, fc = TRUE)
hclr.aov(x, ina)
lr.aov(x, ina)
embed.aov(x, ina)
het.aov(x, ina)
}

\arguments{
\item{x}{
A matrix with the data in Euclidean coordinates, i.e. unit vectors.
}
\item{ina}{
A numerical variable or a factor indicating the group of each vector.
}
\item{fc}{
A boolean that indicates whether a corrected F test should be used or not.
}
}

\details{
The high concentration (hcf.aov), high concentration log-likelihood ratio (hclr.aov),
log-likelihood ratio (lr.aov), embedding approach (embed.aov) or the non equal
concentration parameters approach (het.aov) is used.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
The degree(s) of freedom of the test.
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

Rumcheva P. and Presnell B. (2017). An improved test of equality of mean
directions for the Langevin-von Mises-Fisher distribution.
Australian & New Zealand Journal of Statistics, 59(1): 119--135.

Tsagris M. and Alenazi A. (2022). An investigation of hypothesis testing procedures for circular and spherical mean vectors. Communications in Statistics-Simulation and Computation (Accepted for publication).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
and Giorgos Athineou <gioathineou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hcf.boot}, \link{hcfboot}, \link{hclr.circaov},
}
}

\examples{
x <- rvmf(60, rnorm(3), 15)
ina <- rep(1:3, each = 20)
hcf.aov(x, ina)
hcf.aov(x, ina, fc = FALSE)
lr.aov(x, ina)
embed.aov(x, ina)
het.aov(x, ina)
}


