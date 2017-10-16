\name{Circular-linear correlation}
\alias{circlin.cor}
\title{
Circular-linear correlation
}
\description{
It calculates the squared correlation between a circular and one or more linear variables.
}
\usage{
circlin.cor(theta, x, rads = FALSE)
}
\arguments{
  \item{theta}{
  The circular variable.
}
  \item{x}{
  The linear variable or a matrix containing many linear variables.
}
  \item{rads}{
  If the circualr variable is in rads, this should be TRUE and FALSE otherwise.
}
}

\details{
The squared correlation between a circular and one or more linear variables is calculated.
}
\value{
A matrix with as many rows as linear variables including:
\item{R-squared}{
The value of the squared correlation.
}
\item{p-value}{
The p-value of the zero correlation hypothesis testing.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{circ.cor1}, \link{circ.cor2}, \link{spml.reg}
}
}

\examples{
phi <- rvonmises(50, 2, 20, rads = TRUE)
x <- 2 * phi + rnorm(50)
y <- matrix(rnorm(50 * 5), ncol = 5)
circlin.cor(phi, x, rads = TRUE)
circlin.cor(phi, y, rads = TRUE)
}
\keyword{ Circular-linear correlation }
