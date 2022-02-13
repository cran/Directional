\name{Column-wise uniformity Watson test for circular data}
\alias{colwatsons}
\title{
Column-wise uniformity tests for circular data
}

\description{
Column-wise uniformity tests for circular data.
}

\usage{
colwatsons(u, rads = FALSE)
}

\arguments{
\item{u}{
A numeric matrix containing the circular data which are expressed in radians.
Each column is a different sample.
}
\item{rads}{
A boolean variable. If the data are in radians, put this TRUE.
If the data are expressed in degrees make this FALSE.
}
}

\details{
These tests are used to test the hypothesis that the data come from a circular
uniform distribution.
}

\value{
A matrix with two columns, the value of the test statistic and its associated p-value.
}


\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, pg. 156-157.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{watson}, \link{kuiper}, \link{fishkent}
}
}

\examples{
x <- matrix( rvonmises(n = 50 * 10, m = 2, k = 0), ncol = 10 )
res<-colwatsons(x)
x <- NULL
}



