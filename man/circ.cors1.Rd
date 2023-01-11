\name{Circular correlations between one and many circular variables}
\alias{circ.cors1}
\alias{circ.cors2}

\title{Circular correlations between two circular variables
}

\description{
Circular correlations between two circular variables.
}

\usage{
circ.cors1(theta, phi, rads = FALSE)
circ.cors2(theta, phi, rads = FALSE)
}

\arguments{
\item{theta}{
The first cirular variable expressed in radians, not degrees.
}
\item{phi}{
The other cirular variable. In the case of "circ.cors1" this is a matrix with many
circular variables. In either case, the values must be in radians, not degrees.
}
\item{rads}{
If the data are expressed in rads, then this should be TRUE. If the data are in degrees, then this is FALSE.
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

Jammalamadaka, S. R. and Sarma, Y. R. (1988). A correlation coefficient for angular variables. Statistical Theory and Data Analysis, 2:349--364.

Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
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
circ.cors1(y, x, rads = TRUE)
}

