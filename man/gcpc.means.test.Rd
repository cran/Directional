\name{Two sample location test for circular data under the GCPC distribution}
\alias{gcpc.means.test}

\title{
Two sample location test for circular data under the GCPC distribution.
}

\description{
Two sample location test for circular data under the GCPC distribution.
}

\usage{
gcpc.means.test(u1, u2, rads = FALSE)
}

\arguments{
\item{u1}{
A numeric vector containing circular data for the first sample.
}
\item{u2}{
A numeric vector containing circular data for the second sample.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
}

\details{
The log-likelihood ratio test compares the location parameter of two independent samples,
assuming that both samples are drawn from populations that follow the GCPC distribution.
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

Alzeley O. and Tsagris M. (2026). On the generalized circular projected Cauchy distribution.
https://arxiv.org/abs/2603.04030.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{hcf.circaov}, \link{spcauchy2test}
}
}

\examples{
u1 <- rgcpc(50, omega = 2, g = 5, rho = 0.5, rads = TRUE)
u2 <- rgcpc(50, omega = 2, g = 10, rho = 5, rads = TRUE)
gcpc.means.test(u1, u2, rads = TRUE)
}




