\name{Uniformity test for circular data}
\alias{kuiper}
\alias{watson}

\title{
Uniformity tests for circular data.
}

\description{
Hypothesis tests of uniformity for circular data.
}

\usage{
kuiper(u, rads = FALSE, R = 1)
watson(u, rads = FALSE, R = 1)
}

\arguments{
\item{u}{
A numeric vector containing the circular data, which cna be expressed in degrees or radians.
}
\item{rads}{
A boolean variable. If the data are in radians, put this TRUE. If the data are expressed in degrees make this FALSE.
}
\item{R}{
If R = 1 the asymptotic p-value will be calcualted. If R is greater than 1 the bootstrap p-value is returned.
}
}

\details{
The high concentration (hcf.circaov), log-likelihood ratio (lr.circaov), embedding approach (embed.circaov) or the non equal concentration parameters approach (het.circaov) is used.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
This is usually the degrees of freedom of the test, but here this is "NA" because the asymptotic based p-value is computed in a different way or because bootstrap was employed.
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
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, pg. 153--55 (Kuiper's test) and pg. 156--157 (Watson's test).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rayleigh}, \link{ptest}, \link{vmf.mle}, \link{rvonmises}
}
}

\examples{
x <- rvonmises(n = 40, m = 2, k = 10)
kuiper(x, rads = TRUE)
watson(x, rads = TRUE)
x <- rvonmises(40, m = 2, k = 0)
kuiper(x, rads = TRUE)
watson(x, rads = TRUE)
}

\keyword{ Circular data}
\keyword{ Uniformity test }
\keyword{ Hypothesis testing }



