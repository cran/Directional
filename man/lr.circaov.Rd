\name{lr.circaov}
\alias{lr.circaov}
\title{
Analysis of variance for circular data using the log-likelihood ratio test. 
}
\description{
Analysis of variance for circular data using the the log-likelihood ratio test. 
}
\usage{
lr.circaov(u, ina, rads = F)
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
\value{
A vector including:
\item{test}{
The value of the test statistic.
}
\item{p-value}{
The p-value of the test.
}
\item{kappa}{
The concentration parameter based on all the data. 
}
}
\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{embed.circaov}, \link{conc.test}, \link{hcf.circaov}, \link{het.circaov}
}
\examples{
x <- rvonmises(100, 2.4, 15)
ina <- rep(1:4,each = 25)
lr.circaov(x, ina, rads = TRUE)
}
\keyword{ Circular data}
\keyword{ Analysis of variance }
\keyword{ Hypothesis testing }


