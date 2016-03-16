\name{circ.summary}
\alias{circ.summary}
\title{
Summary statistics for circular data.
}
\description{
It produces a few summary measures for circular data.
}
\usage{
circ.summary(u, rads = FALSE, plot = TRUE)
}
\arguments{
  \item{u}{
  The circular data.
}
  \item{rads}{
If the data are in rads, then this should be TRUE, otherwise FALSE.
}
  \item{plot}{
If you want to see the data plotted on a cicrle make this TRUE.
}
}
\details{
It returns the circular mean, mean resultant length, variance, standard deviation and concentration parameter. So, basically it returns the estimated values of the parameters of the von Mises distribution.
}

\value{
A list including:
\item{mesos}{
The circular mean direction.
}
\item{confint}{
The 95\% confidence interval for the circular mean direction.
}
\item{kappa}{
The concentration parameter.
}
\item{MRL}{
The mean resultant length.
}
\item{circvariance}{
The circular variance.
}
\item{circstd}{
The circular standard deviation.
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
\code{\link{rvonmises}, \link{vm.kde}
}
}
\examples{
x <- rvonmises(50, 2.5, 15, rads = TRUE)
circ.summary(x, rads = TRUE, plot = TRUE)
}
\keyword{ Summary statistics}
\keyword{ Von Mises distribution }
