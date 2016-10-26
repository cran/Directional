\name{Circurlar correlations between two circular variables}
\alias{circ.cor1}
\alias{circ.cor2}

\title{Circurlar correlations between two circular variables
}

\description{
Circurlar correlations between two circular variables.
}

\usage{
circ.cor1(theta, phi, rads = FALSE)

circ.cor2(theta, phi, rads = FALSE)
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
circ.cor1: Correlation for circular variables using the cosinus and sinus formula of Jammaladaka and SenGupta (1988).

circ.cor2: Correlation for circular variables using the cosinus and sinus formula of Mardia and Jupp (2000).
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
Jammalamadaka, R. S. and Sengupta, A. (2001). Topics in circular statistics. World Scientific.

Jammalamadaka, S. R. and Sarma, Y. R. (1988) . A correlation coefficient for angular variables. Statistical Theory and Data Analysis, 2:349--364.

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
\code{ \link{circlin.cor}, \link{circ.cor2}, \link{spml.reg}
}
}

\examples{
y <- runif(50, 0, 2 * pi)
x <- runif(50, 0, 2 * pi)
circ.cor1(x, y, rads = TRUE)
circ.cor2(x, y, rads = TRUE)
}

\keyword{ Circular correlation type I }
\keyword{ Circular correlation type II }
