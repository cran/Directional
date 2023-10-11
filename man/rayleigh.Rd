\name{Rayleigh's test of uniformity}
\alias{rayleigh}
\title{
Rayleigh's test of uniformity
}

\description{
It checkes whether the data are uniformly distributed on the circle or the (hyper-)sphere.
}

\usage{
rayleigh(x, modif = TRUE, B = 999)
}

\arguments{
\item{x}{
A matrix containing the data, unit vectors.
}
\item{modif}{
If modif is TRUE, the modification as suggested by Jupp (2001) is used.
}
\item{B}{
If B is greater than 1, bootstap calibation os performed. If it is equal to 1, classical theory is used.
}
}

\details{
The Rayleigh test of uniformity is not the best, when there are two antipodal mean directions. In this case it will fail. It is good to test whether there is one mean direction or not. To put it differently, it tests whether the concentration parameter of the Fisher distribution is zero or not.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
The degrees of freedom of the test. If bootstrap was employed this is "NA".
}
\item{p.value}{
The p-value of the test.
}
\item{alternative}{
A character with the alternative hypothesis.
}
\item{method}{
A character with the test used.
}
\item{data.name}{
A character vector with two elements.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Jupp, P. E. (2001). Modifications of the rayleigh and bingham tests for uniformity of directions. Journal of Multivariate Analysis, 77(2): 1-20.

Rayleigh, L. (1919). On the problem of random vibrations, and of random flights in one, two, or three dimensions.
The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science, 37(220): 321--347.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{ptest}, \link{kuiper}, \link{iagesag}
}
}

\examples{
x <- rvmf(100, rnorm(5), 1)  ## Fisher distribution with low concentration
rayleigh(x)
}
