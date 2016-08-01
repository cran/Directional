\name{k-NN algorithm using the arc cosinus distance}
\alias{dirknn}
\title{
k-NN algorithm using the arc cosinus distance
}

\description{
It classifies new observations to some known groups via the k-NN algorithm.
}

\usage{
dirknn(x, xnew, k = 5, ina, type = "S", mesos = TRUE)
}

\arguments{
  \item{x}{
  The data, a numeric matrix with unit vectors.
}
  \item{xnew}{
  The new data whose membership is to be predicted, a numeric matrix with unit vectors.
}
  \item{k}{
  The number of nearest neighbours, set to 5 by default.
}
  \item{ina}{
  A variable indicating the groups of the data x.
}
  \item{type}{
  If type is "S", the standard k-NN algorithm is to be used, else "NS" for the non standard one. See below (details) for more information.
}
  \item{mesos}{
  A boolean variable used only in the case of the non standard algorithm (type="NS"). Should the average of the distances be calculated (TRUE) or not (FALSE)? If it is FALSE, the harmonic mean is calculated.
}
}

\details{
The standard algorithm is to keep the k nearest observations and see the groups of these observations. The new observation is allocated to the most frequent seen group. The non standard algorithm is to calculate the classical mean or the harmonic mean of the k nearest observations for each group. The new observation is allocated to the group with the smallest mean distance.
}
\value{
A vector including:
\item{g}{
The predicted group.
}
}

\references{
course webpage of Howard E. Haber.
http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Giorgos Athineou <athineou@csd.uoc.gr>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dirknn.tune}, \link{vmfda.pred}, \link{mix.vmf}
}
}

\examples{
k <- runif(4, 4, 20)
prob <- c(0.2, 0.4, 0.3, 0.1)
mu <- matrix(rnorm(16), ncol = 4)
mu <- mu / sqrt( rowSums(mu^2) )
da <- rmixvmf(200, prob, mu, k)
nu <- sample(1:200, 180)
x <- da$x[nu, ]
ina <- da$id[nu]
xx <- da$x[-nu, ]
id <- da$id[-nu]
a1 <- dirknn(x, xx, k = 5,  ina, type = "S", mesos = TRUE)
a2 <- dirknn(x, xx, k = 5,ina, type = "NS", mesos = TRUE)
a3 <- dirknn(x, xx, k = 5, ina, type = "S", mesos = FALSE)
a4 <- dirknn(x, xx, k = 5, ina, type = "NS", mesos = FALSE)
b <- vmfda.pred(xx, x, ina)
table(id, a1)
table(id, a2)
table(id, a3)
table(id, a4)
}

\keyword{ Directional k-NN algorithm }
\keyword{ Supervised classification }
