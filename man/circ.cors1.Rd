\name{Circurlar correlations between one and many circular variables}
\alias{circ.cors1}

\title{Circurlar correlations between two circular variables
}

\description{
Circurlar correlations between two circular variables.
}

\usage{
circ.cors1(theta, phi)
}

\arguments{
\item{theta}{
The first cirular variable expressed in radians, not degrees.  
}
\item{phi}{
The other cirular variable. In the case of "circ.cors1" this is a matrix with many 
circular variables. In either case, the values must be in radians, not degrees. 
}
}

\details{
Correlation for circular variables using the cosinus and sinus formula of Jammaladaka and SenGupta (1988).
}

\value{
A matrix with two columns, the correlations and the p-values.
}

\references{
Jammalamadaka, R. S. and Sengupta, A. (2001). Topics in circular statistics. World Scientific.

Jammalamadaka, S. R. and Sarma, Y. R. (1988) . A correlation coefficient for angular variables. Statistical Theory and Data Analysis, 2:349--364.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{spml.reg}
}
}

\examples{
y <- runif(50, 0, 2 * pi)
x <- matrix(runif(50 * 10, 0, 2 * pi), ncol = 10)
circ.cors1(y, x)
}

