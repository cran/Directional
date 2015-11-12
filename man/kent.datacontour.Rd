\name{kent.datacontour}
\alias{kent.datacontour}
\title{
Contour plot of the Kent distribution without any data. 
}
\description{
The contour plot of the Kent distribution on the sphere for some data is produced.  
}
\usage{
kent.datacontour(u)
}
\arguments{
  \item{u}{
  A tow column matrix, where the first column is the latitude and the second comlumn is the longitude.
}
}
\details{
MLE of the parameters of the Kent distribution are calculated, then the contour plot is plotted using these estimates and finally the data are also plotted. 
}
\value{
A plot containing the contours of the distribution along with the data. 
}
\references{
Kent John (1982). The Fisher-Bingham distribution on the sphere. Journal of the Royal Statistical Society, Series B, 44(1): 71-80.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{kent.contour}, \link{kent.mle}, \link{vmf.kerncontour}}
}
\examples{
x = rvmf(100, rnorm(3), 10) 
kent.mle(x) 
y = euclid.inv(x) 
kent.datacontour(y)
vmf.kerncontour(y, thumb = "none")
}
\keyword{ Contour plot}
\keyword{ Kent distribution }

