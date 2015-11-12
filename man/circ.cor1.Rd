\name{circ.cor1}
\alias{circ.cor1}
\title{
Correlation for circular variables using the cosinus and sinus formula of Jammaladaka and SenGupta (1988).
}
\description{
Calculates a correlation between two circular variables.
}
\usage{
circ.cor1(theta, phi, rads = F)
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
Jammalamadaka, R. S. and Sengupta, A. (2001). Topics in circular statistics. World Scientific. 
Jammalamadaka, S. R. and Sarma, Y. R. (1988) . A correlation coefficient for angular variables. Statistical Theory and Data Analysis, 2:349--364.}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{circ.cor2}, \link{circlin.cor}}
}
\examples{
y <- runif(50, 0, 2*pi)
x <- runif(50, 0, 2*pi)
circ.cor1(x, y, rads = TRUE)  
}
\keyword{Circular correlation type I}
