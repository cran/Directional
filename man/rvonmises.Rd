\name{Simulation of random values from some circular distributions}
\alias{rvonmises}
\alias{rwrapcauchy}
\alias{rspml}
\alias{rcircbeta}
\alias{rcircpurka}
\alias{rcircexp}
\alias{rcipc}
\alias{rgcpc}
\title{
Simulation of random values from some circular distributions
}

\description{
Simulation of random values from some circular distributions.
}

\usage{
rvonmises(n, m, k, rads = TRUE)
rwrapcauchy(n, m, rho, rads = TRUE)
rspml(n, mu, rads = TRUE)
rcircbeta(n, m, a, b, rads = TRUE)
rcircpurka(n, m, a, rads = TRUE)
rcircexp(n, lambda, rads = TRUE)
rcipc(n, mu = NULL, omega, g, rads = TRUE)
rgcpc(n, mu = NULL, omega, g, rho, rads = TRUE)
}


\arguments{
\item{n}{
The sample size.
}
\item{m}{
The mean angle expressed in radians or degrees.
}
\item{mu}{
The mean vector of the SPML, CIPC and GCPC in \eqn{R^2}. For the CIPC and GCPC, if this argument is not given, 
then the omega and g must be given.
}
\item{omega}{
The location parameter for the CIPC and the GCPC expressed in radians or degrees.
}
\item{k}{
The concentration parameter of the von Mises distribution.
If k is zero the sample will be generated from the uniform distribution over \eqn{(0, 2\pi)}.
}
\item{g}{
The norm of the mean vector for the CIPC and GCPC, if \code{omega} is given instead of \code{mu}.
}
\item{rho}{
For the wrapped Cauchy distribution, this is the \eqn{\rho} parameter.
For the GCPC distribution this is the eigenvalue parameter, or covariance determinant.
}
\item{a}{
The \eqn{\alpha} parameter of the beta distribution.
}
\item{b}{
The \eqn{\beta} parameter of the beta distribution.
}
\item{lambda}{
The \eqn{\lambda} parameter of the circular (wrapped) exponential distribution.
}
\item{rads}{
If the mean angle is expressed in radians, this should be TRUE and FALSE otherwise.
The simulated data will be expressed in radians or degrees depending on what the mean angle is expressed.
}
}

\details{
For the von Mises distribution, the mean direction is transformed to the Euclidean coordinates (i.e. unit vector)
and then the \code{\link{rvmf}} function is employed. It uses a rejection smapling as suggested by Andrew Wood in 1994.
We have mentioned the description of the algorithm as we found it in Dhillon and Sra in 2003.
Finally, the data are transformed to radians or degrees.

For the wrapped Cauchy and wrapped exponential distributions the function generates Cauchy or exponential values
\eqn{x} and then wrapps them around the circle \eqn{x = x(mod 2 \pi)}. For the circular beta the function has some
extra steps (see Zheng Sun's master thesis).

For the CIPC and GCPC distributions, data are generated from the bivariate Cauchy distribution,
normalized to have unit norm and then transformed to angles.
}

\value{
A vector with the simulated data.
}

\references{
Wood A.T.A. (1994). Simulation of the von Mises Fisher distribution.
Communications in Statistics-Simulation and Computation, 23(1): 157-164.

Dhillon I.S. and Sra S. (2003). Modeling data using directional distributions.
Technical Report TR-03-06, Department of Computer Sciences, The University of Texas at Austin.
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.75.4122&rep=rep1&type=pdf

Zheng Sun (2006). Comparing measures of fit for circular distributions. Master thesis, University of Victoria.
https://dspace.library.uvic.ca/bitstream/handle/1828/2698/zhengsun_master_thesis.pdf;sequence=1

Lai M. (1994). Some results in the statistical analysis of directional data.
Master thesis, University of Hong Kong.

Presnell B., Morrison S.P. and Littell R.C. (1998). 
Projected multivariate linear models for directional data.
Journal of the American Statistical Association, 93(443): 1068--1077.

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: 
Obtained through Maximum Likelihood Characterization. 
The Indian Journal of Statistics, Series A, 53(1): 70--83

Jammalamadaka S.R. and Kozubowski T.J. (2003). A new family of circular models: 
The wrapped Laplace distributions. Advances and Applications in Statistics, 3(1): 77--103.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
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

