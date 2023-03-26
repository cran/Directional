\name{Read a file as a Filebacked Big Matrix}
\alias{read.fbm}
\title{
Read a file as a Filebacked Big Matrix
}

\description{
Read a file as a Filebacked Big Matrix.
}

\usage{
read.fbm(file, select)
}

\arguments{
\item{file}{
The File to read.
}
\item{select}{
Indices of columns to read (sorted).
The length of select will be the number of columns of the resulting FBM.
}
}

\details{
The functions read a file as a Filebacked Big Matrix object. For more information see the "bigstatsr" package.
}

\value{
A Filebacked Big Matrix object.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}

\seealso{
\code{ \link{vmf.mle}, \link{kent.mle} }
}

\examples{
\donttest{
if (require("bigreadr") ) {
dataset <- matrix( runif(100 * 50, 2*pi, 100), ncol = 50 )
write.csv(dataset, "dataset.csv")
x <- Directional::read.fbm("dataset.csv", select = 1:50, header = TRUE)
}
}
}
