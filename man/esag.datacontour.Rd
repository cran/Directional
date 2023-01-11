\name{Contour plot of the ESAH distribution for some data}
\alias{esag.datacontour}
\title{
Contour plot of the ESAG distribution for some data
}

\description{
The contour plot of the ESAG distribution on the sphere for some data is produced.
}

\usage{
esag.datacontour(x, lat = 20, long = 20)
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
MLE of the parameters of the ESAG distribution is performed calculated, then the contour plot is plotted using these estimates and finally the data are also plotted.
}

\value{
A plot containing the contours of the distribution along with the data.
}

\references{
Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.
}

\author{
Michail Tsagris and Christos Adam.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Christos Adam \email{pada4m4@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kent.datacontour}, \link{esag.mle}
}
}

\examples{
mu <- colMeans( as.matrix( iris[,1:3] ) )
gam <- c(1,0.5)
x <- resag(100, mu, gam)
esag.mle(x)
y <- euclid.inv(x)
\donttest{
esag.datacontour(y)
esag.datacontour(y, lat = 30, long = 30)

}
}

