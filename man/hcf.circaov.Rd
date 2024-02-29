\name{Anova for circular data}
\alias{hcf.circaov}
\alias{hclr.circaov}
\alias{lr.circaov}
\alias{embed.circaov}
\alias{het.circaov}

\title{
Analysis of variance for circular data
}

\description{
Analysis of variance for circular data.
}

\usage{
hcf.circaov(u, ina, rads = FALSE)
hclr.circaov(u, ina, rads = FALSE)
lr.circaov(u, ina, rads = FALSE)
het.circaov(u, ina, rads = FALSE)
embed.circaov(u, ina, rads = FALSE)
}

\arguments{
\item{u}{
A numeric vector containing the data.
}
\item{ina}{
A numerical or factor variable indicating the group of each value.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
}

\details{
The high concentration (hcf.circaov), high concentration likelihood ratio (hclr.aov),
log-likelihood ratio (lr.circaov), embedding approach (embed.circaov) or the non equal
concentration parameters approach (het.circaov) is used.
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
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics.
Chicester: John Wiley & Sons.

Rumcheva P. and Presnell B. (2017). An improved test of equality of mean directions for
the Langevin-von Mises-Fisher distribution. Australian & New Zealand Journal of Statistics, 59(1): 119--135.

Tsagris M. and Alenazi A. (2024). An investigation of hypothesis testing procedures for circular
and spherical mean vectors. Communications in Statistics-Simulation and Computation, 53(3): 1387--1408.
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
\code{\link{hclr.aov}, \link{hcfcirc.boot}, \link{hcfcircboot}
}
}

\examples{
x <- rvonmises(100, 2.4, 15)
ina <- rep(1:4,each = 25)
hcf.circaov(x, ina, rads = TRUE)
lr.circaov(x, ina, rads = TRUE)
het.circaov(x, ina, rads = TRUE)
embed.circaov(x, ina, rads = TRUE)
hclr.circaov(x, ina, rads = TRUE)
}




