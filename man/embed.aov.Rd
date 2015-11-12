\name{embed.aov}
\alias{embed.aov}
\title{
Analysis of variance for (hyper-)spherical data using the embedding approach. 
}
\description{
Analysis of variance for (hyper-)spherical data using the embedding approach. 
}
\usage{
embed.aov(x, ina)
}
\arguments{
  \item{x}{
  The data in Euclidean coordinates, i.e. unit vectors. 
}
  \item{ina}{
  A numerical variable or a factor indicating the group of each vector. 
}
}
\details{
The embedding approach is used. 
}
\value{
A vector including:
\item{F}{
The F test statistic. 
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
\seealso{
\code{\link{lr.aov}, \link{hcf.aov}, \link{het.aov}, \link{spherconc.test}
}
}
\examples{
x <- rvmf(60, rnorm(3), 15)
ina <- rep(1:3, each = 20)
embed.aov(x, ina)
}
\keyword{ Analysis of variance }
\keyword{ (Hyper-)spherical data }
\keyword{ Embedding approach }
