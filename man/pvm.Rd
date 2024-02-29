\name{Cumulative distribution function of circular distributions}
\alias{pvm}
\alias{pspml}
\alias{pwrapcauchy}
\alias{pcircpurka}
\alias{pcircbeta}
\alias{pcardio}
\alias{pcircexp}
\alias{pcipc}
\alias{pcipc}
\alias{pgcpc}
\alias{pmmvm}
\title{
Cumulative distribution function of circular distributions
}

\description{
Cumulative probability distribution of circular distributions.
}

\usage{
pvm(u, m, k, rads = FALSE)
pspml(u, mu, rads = FALSE)
pwrapcauchy(u, m, rho, rads = FALSE)
pcircpurka(u, m, a, rads = FALSE)
pcircbeta(u, m, a, b, rads = FALSE)
pcardio(u, m, rho, rads = FALSE)
pcircexp(u, lambda, rads = FALSE)
pcipc(u, omega, g, rads = FALSE)
pgcpc(u, omega, g, rho, rads = FALSE)
pmmvm(u, m, k, N, rads = FALSE)
}

\arguments{
\item{u}{
A numerical value, either in radians or in degrees.
}
\item{m}{
The mean direction of the von Mises and the multi-modal von Mises distribution in radians or in degrees.
}
\item{mu}{
The mean vector, a vector with two values for the "pspml".
}
\item{omega}{
The location parameter of the CIPC and GCPC distributions.
}
\item{g}{
The norm of the mean vector for the CIPC and GCPC distributions.
}
\item{k}{
The concentration parameter, \eqn{\kappa}.
}
\item{lambda}{
The \eqn{\lambda} parameter of the circular exponential distribution. This must be positive.
}
\item{a}{
The \eqn{\alpha} parameter of the circular Purkayastha distribution or the \eqn{\alpha} parameter of the
circular Beta distribution.
}
\item{b}{
The \eqn{\beta} parameter of the circular beta distribution.
}
\item{rho}{
The \eqn{\rho} parameter of the Cardioid, wrapped Cauchy and GCPC distributions.
}
\item{N}{
The number of modes to consider in the multi-modal von Mises distribution.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
}

\details{
This value calculates the probability of u being less than some value \eqn{\theta}.
}

\value{
The probability that of u being less than \eqn{\theta}, where u follows a circular distribution.
}

\references{
Arthur Pewsey, Markus Neuhauser, and Graeme D. Ruxton (2013). Circular Statistics in R.

Barnett M. J. and Kingston R. L. (2024). A note on the Hendrickson-Lattman phase probability
distribution and its equivalence to the generalized von Mises distribution.
Journal of Applied Crystallography, 57(2).

Jammalamadaka S. R. and Kozubowski T. J. (2003). A new family of circular models:
The wrapped Laplace distributions. Advances and Applications in Statistics, 3(1): 77--103.

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through
Maximum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70--83

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution.
Communications in Statistics--Theory and Methods, 19(6): 1973--1986.

Paula F. V., Nascimento A. D., Amaral G. J. and Cordeiro G. M. (2021).
Generalized Cardioid distributions for circular data analysis. Stats, 4(3): 634--649.

Zheng Sun (2009). Comparing measures of fit for circular distributions. MSc Thesis, University of Victoria.
file:///C:/Users/mtsag/Downloads/zhengsun_master_thesis.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{group.gof}, \link{dvm}, \link{dcircexp},
\link{purka.mle}, \link{dcircpurka}, \link{dmmvm}
}
}

\examples{
pvm(1, 2, 10, rads = TRUE)
pmmvm(1, 2, 10, 3, rads = TRUE)
pcircexp(c(1, 2), 2, rads = TRUE)
pcircpurka(2, 3, 0.3)
}

