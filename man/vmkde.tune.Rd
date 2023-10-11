\name{Tuning of the bandwidth parameter in the von Mises kernel}
\alias{vmkde.tune}
\title{
Tuning of the bandwidth parameter in the von Mises kernel for circular data
}

\description{
Tuning of the bandwidth parameter in the von Mises kernel for circular data. Cross validation is used.
}

\usage{
vmkde.tune(u, low = 0.1, up = 1, rads = TRUE)
}

\arguments{
\item{u}{
The data, a numerical vector.
}
\item{low}{
The lower value of h to search.
}
\item{up}{
The lower value of h to search.
}
\item{rads}{
If the data are in radians this should be TRUE and FALSE otherwise.
}
}

\details{
Tuning of the bandwidth parameter in the von Mises kernel for circula data via cross validation. The procedure is fast because an optimiser is used.
}

\value{
A vector including two elements:
\item{Optimal h}{
The best H found.
}
\item{cv}{
The value of the maximised pseudo-likelihood.
}
}

\references{
Taylor C.C. (2008). Automatic bandwidth selection for circular density estimation. 
Computational Statistics & Data Analysis, 52(7), 3493--3500.

Wand M.P. and Jones M.C. (1994). Kernel smoothing. CrC Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Giorgos Athineou <gioathineou@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vm.kde}, \link{vmfkde.tune}, \link{vmf.kde}
}
}

\examples{
u <- rvonmises(100, 2.4, 10, rads = TRUE)
vmkde.tune(u)
}
\keyword{ Kernel density estimate }
\keyword{ von Mises distribution }
\keyword{ Tuning the bandwidth }
