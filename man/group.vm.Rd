\name{group.vm}
\alias{group.vm}
\title{
Summary statistics for grouped circular data.
}
\description{
It produces a few summary measures for grouped circular data.
}

\usage{
group.vm(group, fi, rads = FALSE)
}

\arguments{
\item{group}{
A matrix denoting the classes. Each row consists of two numbers, the lower and upper points of each class.
}
\item{fi}{
The frequency of each class of data.
}
\item{rads}{
If the data are in rads, then this should be TRUE, otherwise FALSE.
}

}
\details{
It returns the circular mean, mean resultant length, variance, standard deviation and concentration parameter. So, basically it returns the estimated values of the parameters
of the von Mises distribution. The mena resultant length though is group corrected.
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
Pewsey Arthur, Markus Neuhauser and Graeme D. Ruxton (2013). Circular statistics in R. Oxford University Press.

Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{circ.summary}, \link{rvonmises}, \link{vm.kde}
}
}

\examples{
x <- rvonmises(200, 3, 10)
a <- circ.summary(x, rads = TRUE, plot = FALSE)
group <- seq(min(x) - 0.1, max(x) + 0.1, length = 6)
y <- cut(x, breaks = group, length = 6)
group <- matrix( c( group[1], rep(group[2:5], each = 2), group[6]), ncol = 2, byrow = TRUE)
fi <- as.vector( table(y) )
b <- group.vm(group, fi, rads = TRUE)
a
b
}

\keyword{ Summary statistics}
\keyword{ Von Mises distribution }
\keyword{ Grouped data }
