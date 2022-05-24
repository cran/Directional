\name{Bootstrap 2-sample mean test for circular data}
\alias{hcfcirc.boot}
\alias{lrcirc.boot}
\alias{hclrcirc.boot}
\alias{embedcirc.boot}
\alias{hetcirc.boot}

\title{
Bootstrap 2-sample mean test for circular data
}

\description{
Bootstrap 2-sample mean test for circular data.
}

\usage{
hcfcirc.boot(u1, u2, rads = TRUE, B = 999)
lrcirc.boot(u1, u2, rads = TRUE, B = 999)
hclrcirc.boot(u1, u2, rads = TRUE, B = 999)
embedcirc.boot(u1, u2, rads = TRUE, B = 999)
hetcirc.boot(u1, u2, rads = TRUE, B = 999)
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
The number of bootstraps to perform.
}
}

\details{
The high concentration (hcfcirc.boot), the log-likelihood ratio test (lrcirc.boot), 
high concentration log-likelihood ratio (hclrcirc.boot), embedding approach (embedcirc.boot), 
or the non equal concentration parameters approach (hetcirc.boot) is used.
}

\value{
A vector including two or three numbers, the test statistic value, the bootstrap 
p-value of the test and the common concentration parameter kappa based on all the data.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. 
Chicester: John Wiley & Sons.

Rumcheva P. and Presnell B. (2017). An improved test of equality of mean directions for the 
Langevin-von Mises-Fisher distribution. Australian & New Zealand Journal of Statistics, 59(1), 119-135.

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
hcfcirc.boot(u1, u2)
}



