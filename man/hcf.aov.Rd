\name{Anova for (hyper-)spherical data}
\alias{hcf.aov}
\alias{lr.aov}
\alias{embed.aov}
\alias{het.aov}

\title{
Analysis of variance for (hyper-)spherical data.
}

\description{
Analysis of variance for (hyper-)spherical data.
}

\usage{
hcf.aov(x, ina, fc = TRUE)

lr.aov(x, ina)

embed.aov(x, ina)

het.aov(x, ina)
}

\arguments{
  \item{x}{
  The data in Euclidean coordinates, i.e. unit vectors.
}
  \item{ina}{
  A numerical variable or a factor indicating the group of each vector.
}
  \item{fc}{
  A boolean that indicates whether a corrected F test should be used or not.
}
}

\details{
The high concentration (hcf.aov), log-likelihood ratio (lr.aov), embedding approach (embed.aov) or the non equal concentration parameters approach (het.aov) is used.}

\value{
A vector including:
\item{test}{
The test statistic value.
}
\item{p-value}{
The p-value of the F test.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{lr.aov}, \link{embed.aov}, \link{het.aov}, \link{spherconc.test}
}
}

\examples{
x <- rvmf(60, rnorm(3), 15)
ina <- rep(1:3, each = 20)
hcf.aov(x, ina)
hcf.aov(x, ina, fc = FALSE)
lr.aov(x, ina)
embed.aov(x, ina)
het.aov(x, ina)
}

\keyword{ Analysis of variance }
\keyword{ (Hyper-)spherical data }
\keyword{ high concentration F test }

