\name{Simulation of random values from the ESAG distribution}
\alias{resag}
\alias{rESAGd}
\title{
Simulation of random values from the ESAG distribution
}

\description{
Simulation of random values from the ESAG distribution.
}

\usage{
resag(n, mu, gam)
rESAGd(n, mu, gam)
}

\arguments{
\item{n}{
A number; how many vectors you want to generate.
}
\item{mu}{
The mean vector the ESAG distribution.
}
\item{gam}{
The \eqn{\gamma} parameters of the ESAG distribution. For the rESAGd this may be NULL in case you want to simulate from the IAG in arbitrady dimensions.
}
}

\details{
A random sample from the ESAG distribution is generated. In case the \eqn{\gamma_s} are zero (or null for the rESAGd), the sample is drawn from the Independent Angular Gaussian (or projected normal). The resag() is designed for the sphere, whereas the rESAGd is designed for the sphere and hyper-sphere.
}

\value{
An \eqn{n \times d} matrix with the simulated unit vectors.
}

\references{
Zehao Yu and Xianzheng Huang (2024). A new parameterization for elliptically symmetric angular Gaussian distributions of arbitrary dimension. Electronic Journal of Statististics, 18(1): 301--334.

Paine P.J., Preston S.P., Tsagris M. and Wood A.T.A. (2018). An Elliptically Symmetric Angular
Gaussian Distribution. Statistics and Computing, 28(3):689--697.


Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris and Zehao Yu.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Zehao Yu \email{zehaoy@email.sc.edu}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{esag.mle}, \link{desag}, \link{spml.mle}, \link{acg.mle}, \link{circ.summary}
}
}

\examples{
m <- colMeans( as.matrix( iris[, 1:3] ) )
y <- resag(1000, m, c(1, 0.5) )
esag.mle(y)
}

\keyword{ ESAG distribution }
\keyword{ simulation }
\keyword{ spherical data }
