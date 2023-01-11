\name{Contour plot of the Kent distribution for some data}
\alias{kent.datacontour}
\title{
Contour plot of the Kent distribution for some data
}

\description{
The contour plot of the Kent distribution on the sphere for some data is produced.
}

\usage{
kent.datacontour(x, lat = 20, long = 20)
}

\arguments{
\item{x}{
A two column matrix, where the first column is the latitude and the second comlumn is the longitude. If the matrix has two columns, it is assumed to have unit vectors and in this case it is turned into latitude and longitude.
}
\item{lat}{
A positive number determing the range of degrees to move left and right from the latitude center. See the example to better understand this argument.
}
\item{long}{
A positive number determing the range of degrees to move up and down from the longitude center. See the example to better
understand this argument.
}
}
\details{
MLE of the parameters of the Kent distribution is performed, then the contour plot is plotted using these estimates and finally the data are also plotted.
}

\value{
A plot containing the contours of the distribution along with the data.
}

\references{
Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society, Series B, 44(1): 71-80.
}

\author{
Michail Tsagris and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kent.contour}, \link{kent.mle}, \link{vmf.kerncontour}
}
}

\examples{
x <- rvmf(100, rnorm(3), 10)
kent.mle(x)
y <- euclid.inv(x)
\donttest{
kent.datacontour(y)
kent.datacontour(y, lat = 30, long = 30)
}
}


