\name{rbingham}
\alias{rbingham}
\title{
Simulation from a Bingham distribution using any symmetric 
matrix A}

\description{
It simulates random values from a Bingham distribution with any given symmetric matrix. 
}
\usage{
rbingham(n, A)
}
\arguments{
  \item{n}{
  The sample size.
}
  \item{A}{
  A symmetric matrix.
}
}
\details{
The eigenvalues are fist calcualted and then Chris Fallaize and Theo Kypraio's code (f.rbing) is used. The resulting simulated data anre then right multiplied by the eigenvectors of the matrix A. 
}
\value{
A matrix with the siumlated data.
}
\references{
Kent J.T., Ganeiber A.M. and Mardia K.V. (2013). A new method to simulate the Bingham and related distributions in directional data analysis with applications http://arxiv.org/pdf/1310.8110v1.pdf 

C.J. Fallaize and T. Kypraios (2014). Exact Bayesian Inference for the Bingham Distribution. Statistics and Computing (To appear). 
http://arxiv.org/pdf/1401.2894v1.pdf
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{f.rbing}, \link{rfb}, \link{rvmf}}
}
\examples{
A = cov(iris[,1:3])
x = rbingham(100, A)
}
\keyword{ Bingham distribution }
\keyword{ Simulation of random values }
