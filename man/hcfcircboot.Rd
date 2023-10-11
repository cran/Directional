\name{Bootstrap ANOVA for circular data}
\alias{hcfcircboot}
\alias{hetcircboot}

\title{
Bootstrap ANOVA for circular data
}

\description{
Bootstrap ANOVA for circular data.
}

\usage{
hcfcircboot(u, ina, rads = TRUE, B = 999)
hetcircboot(u, ina, rads = TRUE, B = 999)
}

\arguments{
\item{u}{
A numeric vector containing the data of all groups.
}
\item{ina}{
The grouping variables. A factor or a numerical vector specifying the groups to which each observation belongs to.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
\item{B}{
The number of bootstraps to perform.
}
}

\details{
The high concentration (hcfcircboot), or the non equal concentration
parameters approach (hetcircboot) is used.
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

Rumcheva P. and Presnell B. (2017). An improved test of equality of mean directions for the
Langevin-von Mises-Fisher distribution. Australian & New Zealand Journal of Statistics, 59(1): 119--135.

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
\code{\link{hcf.circaov}, \link{het.aov}
}
}

\examples{
u1 <- rvonmises(20, 2.4, 5)
u2 <- rvonmises(20, 2.4, 10)
hcfcirc.boot(u1, u2)
}



