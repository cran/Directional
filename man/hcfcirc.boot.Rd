\name{Bootstrap 2-sample mean test for circular data}
\alias{hcfcirc.boot}
\alias{lrcirc.boot}
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
The high concentration (hcfcirc.boot), the log-likelihood ratio test (lrcirc.boot) or the non equal concentration parameters approach (hetcirc.boot) is used.
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
The concentration parameter based on all the data. If the het.circaov is used this argument is not returned.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

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



