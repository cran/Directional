\name{kent.contour}
\alias{kent.contour}
\title{
Contour plot of the Kent distribution without any data
}

\description{
The contour plot of the Kent distribution on the sphere is produced. The user can see how the shape and ovalness change as he/she changes the ovlaness parameter.
}

\usage{
kent.contour(k, b)
}

\arguments{
\item{k}{
The concentration parameter.
}
\item{b}{
The ovalness parameter. It has to be less than k/2 in order for the distribution to be unimodal. Otherwise it is bimodal.
}
}

\details{
The goal of this function is for the user to see hwo the Kent distribution looks like.
}

\value{
A plot containing the contours of the distribution.
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
\code{\link{kent.datacontour}, \link{kent.mle}, \link{vmf.contour}, \link{vmf.kerncontour}
}
}

\examples{
par( mfrow = c(1, 2) )
kent.contour(10, 2)
kent.contour(10, 4)
}

\keyword{ Contour plot}
\keyword{ Kent distribution }
