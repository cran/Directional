\name{Permutation based 2-sample mean test for circular data}
\alias{hcfcirc.perm}
\alias{hetcirc.perm}
\alias{lrcirc.perm}
\alias{hclrcirc.perm}
\alias{embedcirc.perm}

\title{
Permutation based 2-sample mean test for circular data
}

\description{
Permutation based 2-sample mean test for circular data.
}

\usage{
hcfcirc.perm(u1, u2, rads = TRUE, B = 999)
hetcirc.perm(u1, u2, rads = TRUE, B = 999)
lrcirc.perm(u1, u2, rads = TRUE, B = 999)
hclrcirc.perm(u1, u2, rads = TRUE, B = 999)
embedcirc.perm(u1, u2, rads = TRUE, B = 999)
}

\arguments{
\item{u1}{
A numeric vector containing the data of the first sample.
}
\item{u2}{
A numeric vector containing the data of the first sample.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
\item{B}{
The number of permutations to perform.
}
}

\details{
The high concentration (hcfcirc.perm), log-likelihood ratio (lrcirc.perm), 
high concentration log-likelihood ratio (hclrcirc.perm), embedding approach (embedcirc.perm) 
or the non equal concentration parameters approach (hetcirc.perm) is used.
}

\value{
A vector including:
\item{test}{
The value of the test statistic.
}
\item{p-value}{
The p-value of the test.
}
\item{kappa}{
The concentration parameter based on all the data. If the het.circaov is used this argument 
is not returned.
}
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Rumcheva P. and Presnell B. (2017). An improved test of equality of mean directions for 
the Langevin-von Mises-Fisher distribution. Australian & New Zealand Journal of Statistics, 59(1): 119--135.

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
\code{\link{hcf.circaov}, \link{het.aov}
}
}

\examples{
u1 <- rvonmises(20, 2.4, 5)
u2 <- rvonmises(20, 2.4, 10)
hcfcirc.perm(u1, u2)
lrcirc.perm(u1, u2)
}



