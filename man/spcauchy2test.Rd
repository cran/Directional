\name{Two sample location test for (hyper-)spherical data}
\alias{spcauchy2test}
\alias{pkbd2test}

\title{
Analysis of variance for (hyper-)spherical data
}

\description{
Analysis of variance for (hyper-)spherical data.
}

\usage{
spcauchy2test(y1, y2, B = 1)
pkbd2test(y1, y2, B = 1)
}

\arguments{
\item{y1}{
A matrix with the data in Euclidean coordinates, i.e. unit vectors.
}
\item{y2}{
A matrix with the data in Euclidean coordinates, i.e. unit vectors.
}
\item{B}{
The number of bootstraps to perform.
}
}

\details{
A log-likelihood ratio based test for the equality of two location parameters, assuming that the data in each group follow the spherical Cauchy of the Poisson-kernel based distribution. Bootstrap is also offered.
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
Kato S. and McCullagh P. (2020). Some properties of a Cauchy family on the sphere derived from the Mobius transformations. Bernoulli, 26(4): 3224--3248.

Golzy M. and Markatou M. (2020). Poisson kernel-based clustering on the sphere: convergence properties,
identifiability, and a method of sampling. Journal of Computational and Graphical Statistics, 29(4): 758--770.

Tsagris M. (2024). Directional data analysis using the spherical Cauchy and the Poisson-kernel based distribution.
https://arxiv.org/pdf/2409.03292
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hcf.boot}, \link{hcfboot}, \link{hclr.circaov}
}
}

\examples{
mu <- rvmf(2, rnorm(5), 3)
y1 <- rspcauchy(60, mu[1, ], 0.4)
y2 <- rspcauchy(30, mu[2, ], 0.8)
spcauchy2test(y1, y2)
}


