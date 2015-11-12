\name{spher.cor}
\alias{spher.cor}
\title{
Spherical-spherical correlation.
}
\description{
Correlation between two spherical variables. 
}
\usage{
spher.cor(x, y)
}
\arguments{
  \item{x}{
  A spherical variable.
}
  \item{y}{
  A spherical variable
}
}
\details{
A very similar to the classical correlation is calcualted. In addition, a hypothesis test for no correlation is performed. Note, that this is a squared correlation actually, so negative values will never be returned. 
}
\value{
A vector including:
\item{R-squared}{
The value of the squared correlation. 
}
\item{p-value}{
The p-value of the no correlation hypothesis testing. 
}
}
\references{
Kanti V. Mardia and Peter E. Jupp. Directional statistics, pg. 254--255. 
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{spher.reg}, \link{vmf}, \link{circ.cor1}, \link{circ.cor2}}
}
\examples{
x <- rvmf(100, rnorm(3), 10)
y <- rvmf(100, rnorm(3), 10)
spher.cor(x, y) 
}
\keyword{ Squared correlation}
\keyword{ Spherical data }
