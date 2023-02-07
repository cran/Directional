\name{Simulation of random values from a Bingham distribution}
\alias{f.rbing}
\title{
Simulating from a Bingham distribution}
\description{
It simulates from a Bingham distribution using the code suggested by Kent et al. (2013).
}

\usage{
f.rbing(n, lam, fast = FALSE)
}

\arguments{
\item{n}{
Sample size.
}
\item{lam}{
Eigenvalues of the diagonal symmetric matrix of the Bingham distribution.
}
\item{fast}{
If you want a fast, efficient simulation set this to TRUE.
}
}

\details{
The user must have calculated the eigenvalues of the diagonal symmetric matrix of the Bingham distribution. 
The function accepts the q-1 eigenvalues only. This means, that the user must have subtracted the 
lowest eigenvalue from the rest and give the non zero ones. The function uses rejection sampling and 
it was written by Chris Fallaize and Theo Kypraios (University of Nottingham) and kindly offered. 
Any questions on the code can be addressed to one of the two aforementioned people. 
It is slightly different than the one Kent et al. (2013) suggests.
}

\value{
A list including:
\item{X}{
The simulated data.
}
\item{avtry}{
The estimate of M in the rejection sampling. The average number of simulated values before a value is accepted.
If the argument fast is set to TRUE this information will not appear.
}
}

\references{
Kent J. T., Ganeiber A. M. & Mardia K. V. (2018). A new unified approach for the simulation of a wide class of directional distributions. 
Journal of Computational and Graphical Statistics, 27(2): 291--301.

Kent J.T., Ganeiber A.M. and Mardia K.V. (2013).
A new method to simulate the Bingham and related distributions
in directional data analysis with applications.
http://arxiv.org/pdf/1310.8110v1.pdf

Fallaize C. J. & Kypraios T. (2016). Exact bayesian inference for the Bingham distribution. 
Statistics and Computing, 26(1): 349--360.
http://arxiv.org/pdf/1401.2894v1.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{rfb}, \link{rvmf}, \link{rbingham}, \link{rkent}, link{rsop}
}
}

\examples{
x <- f.rbing( 100, c(1, 0.6, 0.1) )
x
}

\keyword{ Bingham distribution }
\keyword{ Rejection sampling }
\keyword{ Random values simulation }
