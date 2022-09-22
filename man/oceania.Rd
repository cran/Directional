\name{Maps of the world and the continents}
\alias{asia}
\alias{africa}
\alias{europe}
\alias{north.america}
\alias{oceania}
\alias{south.america}
\alias{worldmap}
\title{
 maps of the world and the continents
}
\description{
It produces maps of the world and the continents.
}
\usage{
asia(title = "Asia", coords = NULL)
africa(title = "Africa", coords = NULL)
europe(title = "Europe", coords = NULL)
north.america(title = "North America", coords = NULL)
oceania(title = "Oceania", coords = NULL)
south.america(title = "South America", coords = NULL)
worldmap(title = "World map", coords = NULL)
}
\arguments{
\item{title}{
A character vector with the title of the map.
}
\item{coords}{
If you want specific points to appear on the plot give the coordinates as a matrix, where the first column
contains the longitude and the second column contains the latitude, in degrees.
}
}
\details{
Maps of the world or the continents are produced. This are experimental functions and plot the countries with specific
colouring at the moment. More functionalities will be added in the future.
}
\value{
A map of the selected continent or the whole world.
}

\author{
Christos Adam.

R implementation and documentation: Christos Adam \email{pada4m4@gmail.com} and Michail Tsagris.
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{sphereplot}
}
}
\examples{
\donttest{
oceania()
}
}
\keyword{ Median direction }
