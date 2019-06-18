\name{Random values simulation from some circular distributions}
\alias{rvonmises}
\alias{rwrapcauchy}
\title{
Random values simulation from some circular distributions
}

\description{
Random values simulation from some circular distributions.
}

\usage{
rvonmises(n, m, k, rads = TRUE)
rwrapcauchy(n, m, rho, rads = TRUE)
}

\arguments{
\item{n}{
The sample size.
}
\item{m}{
The mean angle expressed in radians or degrees.
}
\item{k}{
The concentration parameter of the von Mises distribution. If k is zero the sample will be generated from the uniform distribution over \eqn{(0, 2\pi)}.
}
\item{rho}{
The \eqn{\rho} parameter of the Wrapped Cauchy distribution.
}
\item{rads}{
If the mean angle is expressed in radians, this should be TRUE and FALSE otherwise. The simulated data will be expressed in radians or degrees depending on what the mean angle is expressed.
}
}

\details{
For the von Mises distribution, the mean direction is transformed to the Euclidean coordinates (i.e. unit vector) and then the rvmf function is employed. It uses a rejection smapling as suggested by Andrew Wood in 1994. I have mentioned the description of the algorithm as I found it in Dhillon and Sra in 2003. Finally, the data are transformed to radians or degrees.

For the wrapped Cauchy distribution the function generates Cauchy values \eqn{x} and then wrapps around the circle \eqn{x = x(mod 2 \pi)}. For the circular beta the function has some extra steps (see Zheng Sun's master thesis).
}

\value{
A vector with the simulated data.
}

\references{
Wood, A. T. (1994). Simulation of the von Mises Fisher distribution. Communications in statistics-simulation and computation, 23(1): 157-164.

Dhillon, I. S., & Sra, S. (2003). Modeling data using directional distributions. Technical Report TR-03-06, Department of Computer Sciences, The University of Texas at Austin.
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.75.4122&rep=rep1&type=pdf

Zheng Sun (2006). Comparing measures of fit for circular distributions. Master thesis, University of Victoria.
https://dspace.library.uvic.ca/bitstream/handle/1828/2698/zhengsun_master_thesis.pdf;sequence=1

Lai, M. (1994). Some results in the statistical analysis of directional data. Master thesis, University of Hong Kong.

}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <gioathineou@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{circ.summary}, \link{rvmf}, \link{racg}
}
}

\examples{
x <- rvonmises(100, 2, 25, rads = TRUE)
circ.summary(x, rads = TRUE)
}

\keyword{ Random values simulation }
\keyword{ von Mises-Fisher distribution }

