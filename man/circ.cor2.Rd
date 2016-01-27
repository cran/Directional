\name{circ.cor2}
\alias{circ.cor2}
\title{
Correlation for circular variables using the cosinus and sinus formula of Mardia and Jupp (2000).
}
\description{
Calculates a correlation between two circular variables.
}
\usage{
circ.cor2(theta, phi, rads = F)
}
\arguments{
  \item{theta}{
  The first cirular variable.
}
  \item{phi}{
The other cirular variable.
}
  \item{rads}{
If the data are expressed in rads, then this should be TRUE. If the data are in degrees, then this is FALSE.
}

}
\details{
Calculates a correlation between two circular variables.
}
\value{
A vector including:
\item{rho}{
The value of the correlation coefficient.
}
\item{p-value}{
The p-value of the zero correlation hypothesis testing.
}
}
\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{circ.cor1}, \link{circlin.cor}
}
}
\examples{
y <- runif(50, 0, 2*pi)
x <- runif(50, 0, 2*pi)
circ.cor2(x, y, rads = TRUE)
}
\keyword{ Circular correlation type II }

