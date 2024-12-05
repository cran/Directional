\name{Two sample location test for (hyper-)spherical data}
\alias{spcauchy2test}
\alias{pkbd2test}
\alias{vmf2test}
\alias{sp2}
\alias{pk2}
\alias{vmf2}

\title{
Two sample location test for (hyper-)spherical data
}

\description{
Two sample location test for (hyper-)spherical data.
}

\usage{
spcauchy2test(y1, y2, B = 1)
pkbd2test(y1, y2, B = 1)
vmf2test(y1, y2, B = 1)
sp2(y1, y2, tol = 1e-6)
pk2(y1, y2, tol = 1e-6)
vmf2(y1, y2, tol = 1e-6)
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
\item{tol}{
The tolerance value at which to terminate the iterations.
}
}

\details{
A log-likelihood ratio based test for the equality of two location parameters, assuming that the data in each group follow the spherical Cauchy of the Poisson kernel-based distribution. Bootstrap is also offered.

For the von Mises-Fisher distribution we do the same, but for the mean direction.

The functions sp2() and pk2() estimate the common location of the two groups assuming unequal concentration parameters. These functions are used the compute the log-likelihood under the null hypothesis.
So does the function vmf2(), but the mean direction.
}

\value{
The result of the spcauchy2test(), pkbd2test() and vmf2test() functions is an "htest"class object. Thus it returns a list including:
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
The result of the sp2(), pk2() and vmf2 functions is a list including:
\item{mu}{
The common location parameter, for both samples, under the null hypothesis.
}
\item{rho1}{
The concentration parameter of the first group, asusming a common location parameter.
}
\item{rho2}{
The concentration parameter of the second group, asusming a common location parameter.
}
\item{kappa11}{
The concentration parameter (assuming the von Mises-Fisher distribution) of the first group, asusming a common location parameter.
}
\item{kappa2}{
The concentration parameter (assuming the von Mises-Fisher distribution) of the second group, asusming a common location parameter.
}
\item{loglik}{
The log-likelihood of the whole sample, asusming a common location (or mean direction) parameter.
}
}

\references{
Kato S. and McCullagh P. (2020). Some properties of a Cauchy family on the sphere derived from the Mobius transformations. Bernoulli, 26(4): 3224--3248.

Golzy M. and Markatou M. (2020). Poisson kernel-based clustering on the sphere: convergence properties,
identifiability, and a method of sampling. Journal of Computational and Graphical Statistics, 29(4): 758--770.

Tsagris M., Papastamoulis P. and Kato S. (2024). Directional data analysis using the spherical Cauchy and the Poisson kernel-based distribution. https://arxiv.org/pdf/2409.03292.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{het.boot}, \link{het.aov}
}
}

\examples{
mu <- rvmf(2, rnorm(5), 3)
y1 <- rspcauchy(60, mu[1, ], 0.4)
y2 <- rspcauchy(30, mu[2, ], 0.8)
spcauchy2test(y1, y2)
}


