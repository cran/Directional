\name{Mle of the Matrix Fisher distribution on SO(3)}
\alias{matrixfisher.mle}
\title{
Mle of the Matrix Fisher distribution on SO(3)
}

\description{
It returns the maximum likelihood estimate of the Matrix Fisher parameter F(3x3).
}

\usage{
matrixfisher.mle(X)
}

\arguments{
\item{X}{
An array containing rotation matrices in SO(3).
}
}

%\details{
%Firstly calculate the sample mean \eqn{\bar(X) =  n^(-1) \sum_{i=1}^{n} X_i} for a given random samples \eqn{X_1, ... ..., X_n} from Matrix Fisher distribution. Then determine \eqn{svd(\bar(X))={\hat(\Delta)} {D_g} {\hat(\Gamma)}}, where \eqn{D_g} is diagonal matrix, and \eqn{svd(\bar(F))= {\Delta} {D_g} {\Gamma}}. According to Prentice (1986),
%\eqn{\hat(\Delta)} and \eqn{\hat(\Gamma)} are the maximum likelihood estimates of \eqn{\Delta} and \eqn{\Gamma} respective%ly.
%}

\value{
The components of \eqn{ svd( \bar(X) ) }.
}

\references{
Prentice M. J. (1986). Orientation statistics without parametric assumptions.Journal of the Royal Statistical Society. Series B: Methodological 48(2).
}

\author{
Anamul Sajib & Chris Fallaize.

R implementation and documentation: Anamul Sajib <pmxahsa@nottingham.ac.uk> & Chris Fallaize.
}

\seealso{
\code{\link{rmatrixfisher}
}
}

\examples{
F <- 10^(-1) * matrix( c(85, 11, 41, 78, 39, 60, 43, 64, 48), ncol = 3 )   ### An arbitrary F matrix
X <- rmatrixfisher(3000, F)
matrixfisher.mle(X)
svd(F)
## Bigger sample is needed to obtain better results
## X <- rmatrixfisher(30000, F)
## matrixfisher.mle(X)
## svd(F)

}

\keyword{ Matrix Fisher distribution }
\keyword{ Maximum likelihood estimation }
